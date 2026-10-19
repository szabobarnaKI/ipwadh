* test_ipwadh_pipeline_N2000.do -- run in Stata 19 from the folder containing ipwadhc.ado, ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do,
* ipwadh_eval.ado (v0.2.0) and data/S3_N2000/episodes.dta.   NOT YET RUN.   Needs about 1-2 GB RAM; exact split = 7.5M rows, stcox ~1-2 min each.
* Expected numbers: ref_pipeline.py (independent Python reference, Breslow, robust cluster SE x sqrt(G/(G-1)) = Stata vce(cluster id)):
*   rows 69,645 (69,572 usable), 1,950 clusters, exact split 7,524,605 rows, |C| = 34,545
*   unweighted            coef -0.37968668  SE 0.04176646
*   IPW exact, evalat(end) coef -0.39199867  SE 0.05370549
*   truncated p1/p99: bounds 0.0163167226 / 5.3806707758  coef -0.39045295  SE 0.05012383
*   diagnostics (adherent pieces, person-time weighted):
*     all : pieces 7498991 time 1234393.8812 mean 0.97962741 p1 0.01631672 p50 0.81917161 p99 5.38067078 max 36.99072910 ESS 0.40536887
*     arm0: pieces 3703276 time  609105.6205 mean 0.98224951 p1 0.01749792 p50 0.80913368 p99 5.41416434 max 19.03252976 ESS 0.45736999
*     arm1: pieces 3795715 time  625288.2607 mean 0.97707317 p1 0.01449170 p50 0.82763723 p99 5.34096789 max 36.99072910 ESS 0.36456781
*   bandpoints(100) : 514,912 rows  mid -0.39218589 / 0.05356386   end -0.39210096 / 0.05356511
*   bandpoints(1000): 113,474 rows  mid -0.39372610 / 0.05313442   end -0.39288409 / 0.05297864
* Marginal truth (log mean ratio) about -0.43; the IPW estimate is about 0.7 SE away (z = +0.71): consistency check only.
clear all
set more off
local d "data/S3_N2000"
local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local sp "id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu)"
timer clear
timer on 1

use "`d'/episodes.dta", clear
assert _N == 69645
ipwadhc, `base' `lists' parts replace
di as res "0. weights estimated"

* ---- 1. exact split, evalat(end), truncation and diagnostics by arm ----------------------------------------------------------------------
ipwadh_split, `sp'
assert r(m) == 1
assert _N == 7524605
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) truncate(1 99) replace
scalar lo_s = r(trunc_lo)
scalar hi_s = r(trunc_hi)
local nl_ = r(n_levels)
foreach x in n_used pt mean p1 p50 p99 max ess {
    local all_`x' = r(`x')
}
forvalues i = 1/2 {
    foreach x in n pt mean p1 p50 p99 max ess {
        local a`i'_`x' = r(`x'_`i')
    }
}
assert `nl_' == 2
assert reldif(lo_s, 0.0163167226) < 1e-4 & reldif(hi_s, 5.3806707758) < 1e-4
assert `all_n_used' == 7498991 & reldif(`all_pt', 1234393.8812) < 1e-6
assert reldif(`all_mean', 0.97962741) < 1e-4 & reldif(`all_p50', 0.81917161) < 1e-4 & reldif(`all_max', 36.99072910) < 1e-4
assert reldif(`all_ess', 0.40536887) < 1e-4
assert `a1_n' == 3703276 & reldif(`a1_pt', 609105.6205) < 1e-6 & reldif(`a1_mean', 0.98224951) < 1e-4
assert reldif(`a1_p99', 5.41416434) < 1e-4 & reldif(`a1_max', 19.03252976) < 1e-4 & reldif(`a1_ess', 0.45736999) < 1e-4
assert `a2_n' == 3795715 & reldif(`a2_pt', 625288.2607) < 1e-6 & reldif(`a2_mean', 0.97707317) < 1e-4
assert reldif(`a2_p99', 5.34096789) < 1e-4 & reldif(`a2_max', 36.99072910) < 1e-4 & reldif(`a2_ess', 0.36456781) < 1e-4
di as res "1. exact split, truncation bounds and diagnostics match Python: ok"

* ---- 2. Cox fits: unweighted, IPW, truncated IPW ----------------------------------------------------------------------------------------
stset t1 if adh == 1, time0(t0) failure(inf)
stcox i.arm, breslow vce(cluster id)
di "unweighted:  coef " %12.8f _b[1.arm] "  SE " %12.8f _se[1.arm]
assert abs(_b[1.arm] - (-0.37968668)) < 1e-4 & abs(_se[1.arm] - 0.04176646) < 1e-4

stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
stcox i.arm, breslow vce(cluster id)
di "IPW exact:   coef " %12.8f _b[1.arm] "  SE " %12.8f _se[1.arm]
assert abs(_b[1.arm] - (-0.39199867)) < 2e-4 & abs(_se[1.arm] - 0.05370549) < 2e-4
local z = (_b[1.arm] + 0.43) / _se[1.arm]
di as txt "  z versus marginal truth -0.43: " %5.2f `z'

stset t1 [pw=_ipw_pswt] if adh == 1, time0(t0) failure(inf)
stcox i.arm, breslow vce(cluster id)
di "truncated:   coef " %12.8f _b[1.arm] "  SE " %12.8f _se[1.arm]
assert abs(_b[1.arm] - (-0.39045295)) < 2e-4 & abs(_se[1.arm] - 0.05012383) < 2e-4
di as res "2. Cox fits (unweighted, IPW, truncated) match Python: ok"

* ---- 3. banded splits, evalat mid and end ------------------------------------------------------------------------------------------------
foreach m in 100 1000 {
    use "`d'/episodes.dta", clear
    quietly ipwadhc, `base' `lists' parts replace
    ipwadh_split, `sp' bandpoints(`m')
    if `m' == 100  local nrow = 514912
    if `m' == 1000 local nrow = 113474
    assert _N == `nrow'
    foreach ev in mid end {
        ipwadh_eval, time0(t0) time1(t1) evalat(`ev') adherence(adh) replace
        stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
        stcox i.arm, breslow vce(cluster id)
        di "bandpoints(`m') `ev':  coef " %12.8f _b[1.arm] "  SE " %12.8f _se[1.arm]
        if `m' == 100 & "`ev'" == "mid"  assert abs(_b[1.arm] - (-0.39218589)) < 2e-4 & abs(_se[1.arm] - 0.05356386) < 2e-4
        if `m' == 100 & "`ev'" == "end"  assert abs(_b[1.arm] - (-0.39210096)) < 2e-4 & abs(_se[1.arm] - 0.05356511) < 2e-4
        if `m' == 1000 & "`ev'" == "mid" assert abs(_b[1.arm] - (-0.39372610)) < 2e-4 & abs(_se[1.arm] - 0.05313442) < 2e-4
        if `m' == 1000 & "`ev'" == "end" assert abs(_b[1.arm] - (-0.39288409)) < 2e-4 & abs(_se[1.arm] - 0.05297864) < 2e-4
    }
}
di as res "3. banded splits (m = 100, 1000; mid, end) match Python: ok"
timer off 1
timer list 1
di as res "test_ipwadh_pipeline_N2000.do: all assertions passed"
