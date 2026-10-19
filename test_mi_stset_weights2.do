* test_mi_stset_weights2.do  --  round 2 (NOT YET RUN).  Round 1 (output.log, 2026-10-11) showed:
*   - mi stset REFUSES an imputed or a passive weight (rc 459), wide and flong.
*   - an UNREGISTERED copy is accepted by mi stset; re-registering it (passive) afterwards makes  mi estimate: stcox  fail with 459 (the stset basis is re-checked).
*   - mi unregister of a passive variable: values survive in flong, are reset to the m=0 values in wide.
*   - C5 of round 1 was a bug of the test (registered regular first, then register imputed -> rc 110; m=0 missing propagated to all m).
* Round 2 asks:
*   C6a flong: UNREGISTERED ("super-varying") weight, m=0 value 1 (non-missing), m>0 values exp(0.2 x); fill first, then mi stset -> does mi stset work, and does mi estimate use the m-specific weights?
*   C6b flong: same, but mi stset FIRST (weight 1 everywhere), the m>0 values are filled afterwards.
*   C7  wide : the same as C6a (expected: the m>0 values are reset, i.e. impossible in wide).
*   C5f fixed workaround c: unregistered placeholder -> mi stset -> m=0 weight missing -> mi register imputed -> m>0 filled -> mi estimate (expected: 459 like C4).
*   C8  NO stset: weighted mean with an IMPUTED weight:  mi estimate: mean x [pw=w_imp] if arm==1  vs the mean over m of the per-m weighted means (the baseline-table route).
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
    gen double x = rnormal(0, 1) + 0.3 * arm
    replace x = . if runiform() < 0.25
    drop T time d
    sort id row
    mi set `style'
    mi register imputed x
    mi register regular id arm row t0 t1 fail
    quietly mi impute regress x = arm t1 fail, add(5) rseed(20261011)
end

capture program drop recap
program define recap
    args text
    global ncase = $ncase + 1
    global r$ncase `"`text'"'
end

capture program drop try_stset
program define try_stset, rclass
    args wvar
    capture noisily mi stset t1 [pw=`wvar'], time0(t0) failure(fail)
    return scalar rc = _rc
end

capture program drop est_check
program define est_check
    args label wref
    capture noisily mi estimate, post: stcox i.arm, breslow vce(cluster id)
    local rc = _rc
    if `rc' {
        recap "`label': mi stset OK, but mi estimate: stcox FAILED, rc = `rc'"
        exit
    }
    local bmi = _b[1.arm]
    local sum = 0
    local bmin = .
    local bmax = .
    forvalues m = 1/5 {
        preserve
        quietly mi extract `m', clear
        quietly stset t1 [pw=`wref'], time0(t0) failure(fail)
        quietly stcox i.arm, breslow vce(cluster id)
        local b = _b[1.arm]
        local sum = `sum' + `b'
        if missing(`bmin') | `b' < `bmin' local bmin = `b'
        if missing(`bmax') | `b' > `bmax' local bmax = `b'
        restore
    }
    local bref = `sum' / 5
    local dif = abs(`bmi' - `bref')
    local spread = `bmax' - `bmin'
    local verdict "MISMATCH (mi estimate does NOT use the m-specific weights)"
    if `dif' < 1e-8 & `spread' > 1e-6 local verdict "OK: mi estimate uses the m-specific weights"
    if `dif' < 1e-8 & `spread' <= 1e-6 local verdict "agrees, but weights hardly differ across m (uninformative)"
    local line = "`label': mi stset OK;  mi estimate log HR " + string(`bmi', "%10.6f") + "  reference mean " + string(`bref', "%10.6f") + "  diff " + string(`dif', "%9.2e") + "  spread over m " + string(`spread', "%9.2e") + "  -> `verdict'"
    recap `"`line'"'
end

capture program drop show_m
program define show_m
    args wvar
    foreach m in 0 1 2 {
        quietly mi xeq `m': summarize `wvar'
        di as txt "    `wvar' in m=`m': N = " r(N) "  mean = " %9.5f r(mean) "  sd = " %9.5f r(sd)
    }
end

* ===== C6a / C6b / C7 / C5f ====================================================================================================
foreach style in flong wide {

    di as res _n "################  style `style'  ################"

    * ----- C6a / C7: unregistered varying weight, fill first, then mi stset ---------------------------------------------
    local tag = cond("`style'" == "flong", "C6a", "C7 ")
    di as res _n "--- `style' `tag': UNREGISTERED weight, m=0 value 1, m>0 = exp(0.2 x); fill first, then mi stset"
    mkdata `style'
    gen double w_sv = 1
    mi xeq 1/5: replace w_sv = exp(0.2 * x)
    capture noisily mi update
    di as txt "  after mi update (are the m>0 values kept?):"
    show_m w_sv
    capture noisily mi varying
    try_stset w_sv
    local rc_set = r(rc)
    if `rc_set' {
        recap "`style' `tag' unregistered varying weight, fill then stset: mi stset REFUSED, rc = `rc_set'"
    }
    else {
        est_check "`style' `tag' unregistered varying weight, fill then stset" w_sv
    }

    * ----- C6b: mi stset first, fill afterwards (flong only) -------------------------------------------------------------
    if "`style'" == "flong" {
        di as res _n "--- `style' C6b: UNREGISTERED weight, mi stset first (weight 1), m>0 filled afterwards"
        mkdata `style'
        gen double w_sv = 1
        try_stset w_sv
        local rc_set = r(rc)
        if `rc_set' {
            recap "`style' C6b unregistered weight, stset then fill: mi stset REFUSED, rc = `rc_set'"
        }
        else {
            mi xeq 1/5: replace w_sv = exp(0.2 * x)
            capture noisily mi update
            di as txt "  after the fill and mi update:"
            show_m w_sv
            est_check "`style' C6b unregistered weight, stset then fill" w_sv
        }
    }

    * ----- C5f: workaround c fixed -----------------------------------------------------------------------------------------
    di as res _n "--- `style' C5f: unregistered placeholder -> mi stset -> m=0 weight missing -> mi register imputed -> m>0 filled"
    mkdata `style'
    gen double w_ph = 1
    try_stset w_ph
    local rc_set = r(rc)
    if `rc_set' {
        recap "`style' C5f workaround c: mi stset REFUSED for the unregistered placeholder, rc = `rc_set'  <-- unexpected"
    }
    else {
        quietly mi xeq 0: replace w_ph = .
        capture noisily mi register imputed w_ph
        local rc_reg = _rc
        quietly mi xeq 1/5: replace w_ph = exp(0.2 * x)
        capture noisily mi update
        di as txt "  after the fill and mi update:"
        show_m w_ph
        recap "`style' C5f workaround c: stset rc 0, register imputed rc `rc_reg'"
        est_check "`style' C5f workaround c (mi estimate)" w_ph
    }
}

* ===== C8: weighted mean with an imputed weight, no stset ====================================================================
foreach style in wide flong {
    di as res _n "--- `style' C8: mi estimate: mean x [pw=w_imp] if arm==1  (imputed weight, m=0 missing, no stset)"
    mkdata `style'
    gen double w_imp = .
    mi register imputed w_imp
    mi xeq 1/5: replace w_imp = exp(0.2 * x)
    capture noisily mi estimate, post: mean x [pw=w_imp] if arm == 1
    local rc = _rc
    if `rc' {
        recap "`style' C8  mi estimate: mean x [pw=w_imp]: FAILED, rc = `rc'"
    }
    else {
        local bmi = _b[x]
        local sum = 0
        local bmin = .
        local bmax = .
        forvalues m = 1/5 {
            preserve
            quietly mi extract `m', clear
            quietly mean x [pw=w_imp] if arm == 1
            local b = _b[x]
            local sum = `sum' + `b'
            if missing(`bmin') | `b' < `bmin' local bmin = `b'
            if missing(`bmax') | `b' > `bmax' local bmax = `b'
            restore
        }
        local bref = `sum' / 5
        local dif = abs(`bmi' - `bref')
        local spread = `bmax' - `bmin'
        local verdict "MISMATCH"
        if `dif' < 1e-8 local verdict "OK: mi estimate: mean uses the m-specific weights"
        local line = "`style' C8  mi estimate: mean x [pw=w_imp] if arm==1: " + string(`bmi', "%10.6f") + "  manual mean over m " + string(`bref', "%10.6f") + "  diff " + string(`dif', "%9.2e") + "  spread over m " + string(`spread', "%9.2e") + "  -> `verdict'"
        recap `"`line'"'
    }
}

di as res _n "=================================== RECAP ==================================="
forvalues i = 1/$ncase {
    di as txt "${r`i'}"
}
di as res _n "test_mi_stset_weights2.do finished: send me the WHOLE log."
