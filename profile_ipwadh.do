* profile_ipwadh.do -- WHERE DOES THE TIME GO?  Run in Stata 19 from the folder with ipwadhc.ado (0.3.0), ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do,
* ipwadh_eval.ado and data/S3_N2000/episodes.dta.     Runtime about 5-8 minutes.  Send me the whole log.
* No ado file is changed: the do-file loads ipwadhc, then REPLACES its subprogram _ipwadhc_fit by an instrumented copy (same code + timers + one printed line per Poisson fit).
* Question: is Stata's -poisson- the rate-limiting step of ipwadhc / of the whole pipeline?  If it is less than about 10-15% of the total, a home-made Mata Poisson fit is not worth it.
clear all
set more off
discard
local d "data/S3_N2000"
local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local sp "id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu)"

* load the ado file (error 100 = syntax, but the program and its subprograms are now in memory), then replace the fit subprogram
capture ipwadhc
capture program drop _ipwadhc_fit
program define _ipwadhc_fit
    args dep rhs touse len strat nq newvar what popts
    quietly gen double `newvar' = .
    quietly levelsof `strat' if `touse', local(lv)
    tempvar mu
    foreach l of local lv {
        quietly count if `touse' & `strat' == `l' & `dep' == 1
        local nev = r(N)
        if `nev' == 0 {
            di as err "no events for the model of `what' in stratum `l'"
            exit 2000
        }
        quietly count if `touse' & `strat' == `l'
        local nstr = r(N)
        local po `"`popts'"'
        local usedtight 0
        if `"`po'"' == "" {
            local po "tolerance(1e-12) ltolerance(1e-14) nrtolerance(1e-12)"
            local usedtight 1
        }
        timer clear 96
        timer on 96
        capture `nq' poisson `dep' `rhs' if `touse' & `strat' == `l', exposure(`len') `po'
        local fallback 0
        if `usedtight' & (_rc | (_rc == 0 & e(converged) == 0)) {
            local fallback 1
            capture `nq' poisson `dep' `rhs' if `touse' & `strat' == `l', exposure(`len')
        }
        timer off 96
        if _rc {
            local rc = _rc
            di as err "Poisson model for `what' failed in stratum `l' (r(`rc'))"
            exit `rc'
        }
        local ic = e(ic)
        local npar = colsof(e(b))
        timer clear 95
        timer on 95
        quietly predict double `mu' if `touse' & `strat' == `l', n
        quietly replace `newvar' = `mu' / `len' if `touse' & `strat' == `l'
        quietly drop `mu'
        timer off 95
        quietly timer list 96
        local t96 = r(t96)
        quietly timer list 95
        local t95 = r(t95)
        global ipw_tpois = $ipw_tpois + `t96'
        global ipw_tpred = $ipw_tpred + `t95'
        global ipw_nfit  = $ipw_nfit + 1
        di as txt %-52s abbrev("`what'", 52) " str `l'  N " %7.0fc `nstr' "  events " %6.0fc `nev' "  par " %3.0f `npar' "  iter " %3.0f `ic' ///
           "  poisson " as res %7.3f `t96' as txt " s  predict " as res %6.3f `t95' as txt " s" cond(`fallback', "  (FALLBACK to default tolerances)", "")
    }
end

capture program drop prof_run
program define prof_run
    * prof_run label : poissonopts-string   (runs ipwadhc with parts and prints the Poisson share)
    args label popts
    global ipw_tpois = 0
    global ipw_tpred = 0
    global ipw_nfit = 0
    timer clear 91
    timer on 91
    if `"`popts'"' == "" ipwadhc, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu) stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker) parts replace
    else ipwadhc, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu) stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker) parts replace poissonopts(`popts')
    timer off 91
    quietly timer list 91
    local tt = r(t91)
    di as res _n "[`label']  ipwadhc total " %7.2f `tt' " s;  " $ipw_nfit " Poisson fits take " %7.2f $ipw_tpois " s (" %4.1f 100*$ipw_tpois/`tt' "%);  predict/replace " %6.2f $ipw_tpred " s"
    global ipw_lastt = `tt'
    global ipw_lastp = $ipw_tpois
end

* ===== 1. the full pipeline at N = 2000, phase by phase ========================================================================================
use "`d'/episodes.dta", clear
di as res _n "=== 1. pipeline, tight tolerances (the default) ==="
timer clear
prof_run "tight (default)" ""
scalar t_c = $ipw_lastt
scalar t_p = $ipw_lastp
timer on 2
quietly ipwadh_split, `sp'
timer off 2
local nsplit = _N
timer on 3
quietly ipwadh_eval, time0(t0) time1(t1) adherence(adh) truncate(1 99) replace
timer off 3
timer on 4
quietly stset t1 [pw=_ipw_pswt] if adh == 1, time0(t0) failure(inf)
quietly stcox i.arm, breslow vce(cluster id)
timer off 4
quietly timer list 2
scalar t_s = r(t2)
quietly timer list 3
scalar t_e = r(t3)
quietly timer list 4
scalar t_x = r(t4)
scalar t_tot = t_c + t_s + t_e + t_x
di as res _n "PHASES (rows after split: " %12.0fc `nsplit' ")"
di as txt "  ipwadhc total               " %9.2f t_c " s  " %5.1f 100*t_c/t_tot "%"
di as txt "    of which Poisson fits     " %9.2f t_p " s  " %5.1f 100*t_p/t_tot "% of the pipeline"
di as txt "  ipwadh_split                " %9.2f t_s " s  " %5.1f 100*t_s/t_tot "%"
di as txt "  ipwadh_eval                 " %9.2f t_e " s  " %5.1f 100*t_e/t_tot "%"
di as txt "  stset + stcox (one model)   " %9.2f t_x " s  " %5.1f 100*t_x/t_tot "%"
di as txt "  pipeline total              " %9.2f t_tot " s"

* ===== 2. cost of the tight tolerances: same run with Stata's default convergence criteria ================================================
use "`d'/episodes.dta", clear
di as res _n "=== 2. Stata's default tolerances (poissonopts(iterate(100))) ==="
prof_run "default tolerances" "iterate(100)"
di as txt "  tight (default): Poisson " %7.2f t_p " s;   default tolerances: Poisson " %7.2f $ipw_lastp " s"

* ===== 3. scaling: 5 x the persons (ids shifted), ipwadhc only ======================================================================================
di as res _n "=== 3. scaling of ipwadhc with the number of rows (copies of the data with new ids) ==="
tempfile one
use "`d'/episodes.dta", clear
quietly save `one', replace
foreach k in 1 5 10 {
    use "`d'/episodes.dta", clear
    if `k' > 1 {
        forvalues j = 2/`k' {
            quietly append using `one'
            quietly replace id = id + `=(`j'-1)*100000' if _n > (`j'-1) * 69645
        }
    }
    quietly count
    local nr = r(N)
    global ipw_tpois = 0
    global ipw_tpred = 0
    global ipw_nfit = 0
    timer clear 91
    timer on 91
    quietly ipwadhc, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu) stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker) parts replace
    timer off 91
    quietly timer list 91
    di as txt "  rows " %9.0fc `nr' ":  ipwadhc " %8.2f r(t91) " s,  Poisson fits " %8.2f $ipw_tpois " s (" %4.1f 100*$ipw_tpois/r(t91) "%),  " $ipw_nfit " fits"
}

discard
di as res _n "profile_ipwadh.do finished. Send me the whole log."
