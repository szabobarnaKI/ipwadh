* test_mi_build_smoke.do  --  smoke test of the ASSEMBLY of the in-memory mi split data (NOT YET RUN; written without access to Stata).  Stata 16+, synthetic data, no files.
* What it certifies (all pieces of PLAN_mi_splitdata.md section 4.1 that are not ipwadh code):
*   R2  after  mi extract m, clear  the rows come in the SAME order as in m=0 (id row t0 t1 arm identical) for wide, flong and mlong source data; is _mi_id kept by mi extract?
*   B1  a wide mi dataset can be BUILT FROM SCRATCH in a frame from a plain pieces dataset:  mi set wide -> mi register imputed -> mi set M = 5 -> st_store into _k_<var> variables
*       (per-m imputed covariate x and a per-m "person-time weight" _ipw_ptw = exp(0.2 x) * piece length, looked up by the row index _ipw_orig of the unsplit data);
*   B2  the m>0 values survive  mi update  (m=0 of _ipw_ptw is missing);
*   B3  the person-time weighted baseline table pieces:  mi estimate: mean x [pw=_ipw_ptw], over(arm)   and   mi estimate: proportion g [pw=_ipw_ptw], over(arm)
*       equal the mean over m of the same commands run on independently rebuilt per-m data (max abs difference of the coefficient vectors < 1e-10).
* "Pieces" are emulated: every unsplit row is cut into two equal halves (expand 2), which is all the lookup needs (_ipw_orig).
clear all
set more off
version 16
global ncase = 0

capture program drop mkdata
program define mkdata
    args style
    clear
    set seed 4711
    set obs 200
    gen long id = _n
    gen byte arm = (runiform() < .5)
    gen double T = rexponential(8) * exp(-0.4 * arm)
    gen double time = min(T, 10)
    gen byte d = (T < 10)
    expand 2
    bysort id: gen byte row = _n
    gen double t0 = cond(row == 1, 0, time/2)
    gen double t1 = cond(row == 1, time/2, time)
    gen byte fail = (row == 2 & d == 1)
    gen byte g = 1 + (t1 > 2) + (t1 > 5)
    gen double x = rnormal(0, 1) + 0.3 * arm
    replace x = . if runiform() < 0.25
    drop T time d
    sort id row
    mi set `style'
    mi register imputed x
    mi register regular id arm row t0 t1 fail g
    quietly mi impute regress x = arm t1 fail, add(5) rseed(20261011)
end

capture program drop recap
program define recap
    args text
    global ncase = $ncase + 1
    global r$ncase `"`text'"'
end

foreach style in wide flong mlong {

    di as res _n "################  source style `style'  ################"
    mkdata `style'
    quietly mi describe
    local M = r(M)
    capture frame drop F0
    frame copy default F0

    * ----- R2: row order and _mi_id after mi extract ------------------------------------------------------------------------------
    frame F0 {
        quietly mi extract 0, clear
        quietly count
        local N0 = r(N)
        capture confirm variable _mi_id
        local hasid0 = (_rc == 0)
        mata: ID0 = st_data(., "id row t0 t1 arm")
    }
    mata: X = J(`N0', `M', .)
    local allsame 1
    forvalues m = 1/`M' {
        capture frame drop Fm
        frame copy default Fm
        frame Fm {
            quietly mi extract `m', clear
            mata: st_local("same", strofreal(rows(st_data(., "id")) == rows(ID0) ? all(st_data(., "id row t0 t1 arm") :== ID0) : 0))
            mata: X[., `m'] = st_data(., "x")
            mata: st_local("nmiss", strofreal(missing(X[., `m'])))
        }
        frame drop Fm
        if `same' != 1 local allsame 0
        di as txt "  m=`m': same row order/structure as m=0: `same';  missing x: `nmiss'"
    }
    recap "`style' R2: mi extract keeps the row order of m=0 in all m: `allsame';  _mi_id kept by mi extract 0: `hasid0'"

    * ----- B1: build the pieces frame and the wide mi dataset from scratch ----------------------------------------------------------
    capture frame drop P
    frame copy F0 P
    frame P {
        capture drop _mi_id
        gen long _ipw_orig = _n
        quietly expand 2
        bysort _ipw_orig: gen byte _ipw_piece = _n
        sort _ipw_orig _ipw_piece
        gen double len = (t1 - t0) / 2
        gen double _ipw_ptw = .
        mi set wide
        mi register imputed x _ipw_ptw
        mi register regular id arm row t0 t1 fail g len _ipw_orig _ipw_piece
        mi set M = `M'
        forvalues m = 1/`M' {
            mata: st_store(., "_`m'_x", X[st_data(., "_ipw_orig"), `m'])
            mata: st_store(., "_`m'__ipw_ptw", exp(0.2 :* X[st_data(., "_ipw_orig"), `m']) :* st_data(., "len"))
        }
        capture noisily mi update
        di as txt "  B2: after mi update, _ipw_ptw in m=0 / 1 / 2:"
        foreach m in 0 1 2 {
            quietly mi xeq `m': summarize _ipw_ptw
            di as txt "    m=`m': N = " r(N) "  mean = " %9.5f r(mean)
        }
        capture noisily mi describe

        * ----- B3: person-time weighted baseline-table pieces, by arm -------------------------------------------------------------------
        capture noisily mi estimate, post: mean x [pw=_ipw_ptw], over(arm)
        local rc_mean = _rc
        matrix bmi_mean = e(b)
        capture noisily mi estimate, post: proportion g [pw=_ipw_ptw], over(arm)
        local rc_prop = _rc
        matrix bmi_prop = e(b)
    }

    * reference: rebuild every m independently from the ORIGINAL mi data (no X matrix, no _k_ variables)
    matrix bref_mean = J(1, colsof(bmi_mean), 0)
    matrix bref_prop = J(1, colsof(bmi_prop), 0)
    forvalues m = 1/`M' {
        capture frame drop Fm
        frame copy default Fm
        frame Fm {
            quietly mi extract `m', clear
            quietly expand 2
            gen double w = exp(0.2 * x) * (t1 - t0) / 2
            quietly mean x [pw=w], over(arm)
            matrix bm = e(b)
            quietly proportion g [pw=w], over(arm)
            matrix bp = e(b)
        }
        frame drop Fm
        matrix bref_mean = bref_mean + bm / `M'
        matrix bref_prop = bref_prop + bp / `M'
    }
    mata: st_numscalar("dmean", max(abs(st_matrix("bmi_mean") - st_matrix("bref_mean"))))
    mata: st_numscalar("dprop", max(abs(st_matrix("bmi_prop") - st_matrix("bref_prop"))))
    local v1 = cond(`rc_mean' == 0 & scalar(dmean) < 1e-10, "OK", "FAIL")
    local v2 = cond(`rc_prop' == 0 & scalar(dprop) < 1e-10, "OK", "FAIL")
    local line = "`style' B1/B3: wide mi dataset built from scratch;  mi estimate: mean, over(arm): rc `rc_mean', max diff to per-m reference " + string(scalar(dmean), "%9.2e") + " -> `v1';  proportion: rc `rc_prop', max diff " + string(scalar(dprop), "%9.2e") + " -> `v2'"
    recap `"`line'"'
    capture frame drop P
    capture frame drop F0
}

di as res _n "=================================== RECAP ==================================="
forvalues i = 1/$ncase {
    di as txt "${r`i'}"
}
di as res _n "test_mi_build_smoke.do finished: send me the WHOLE log."
