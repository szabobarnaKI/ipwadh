* test_ipwadh_eval_trunc.do -- run in Stata 19 from the folder containing ipwadhc.ado, ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do,
* ipwadh_eval.ado (v0.2.0) and data/.   NOT YET RUN.   Data: S3_N200, exact split (77,718 pieces), evalat(end).
* Expected numbers: ref_piece_weights.py section (6) (independent Python reference; same weights up to the Poisson optimizer tolerance ~6e-7):
*   person-time weighted p1/p99 of the combined weight over adherent pieces: 0.0092329060 / 6.0235830652;  Cox on truncated weights: coef -0.4586770716, SE 0.1809860006
*   all  : pieces 75182  time 118443.1655  mean 0.97757218  p1 0.00923291 p50 0.81501708 p99 6.02358307 max 15.82371047 ESS 0.45611578
*   arm 0: pieces 37043  time  58319.5651  mean 1.12640776  p1 0.00634033 p50 0.88543061 p99 6.54839498 max 14.06583325 ESS 0.47700588
*   arm 1: pieces 38139  time  60123.6005  mean 0.83320247  p1 0.02272328 p50 0.76110270 p99 4.57495231 max 15.82371047 ESS 0.44863962
clear all
set more off
local d "data/S3_N200"
local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local sp "id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu)"

use "`d'/episodes.dta", clear
ipwadhc, `base' `lists' parts replace
ipwadh_split, `sp'
local m_ = r(m)
assert `m_' == 1
assert _N == 77718

* ---- 1. run with truncation and diagnostics by arm; copy ALL returned values into locals at once ----------------------------------
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) truncate(1 99) replace
local wv_    = "`r(weightvar)'"
scalar lo_s = r(trunc_lo)             // scalars: full double precision (a local macro keeps ~16 digits)
scalar hi_s = r(trunc_hi)
local slo_   = r(share_lo)
local shi_   = r(share_hi)
local nl_    = r(n_levels)
foreach x in n_used pt mean min p1 p50 p99 max ess {
    local all_`x' = r(`x')
}
forvalues i = 1/2 {
    local lev`i' = r(level_`i')
    foreach x in n pt mean p1 p50 p99 max ess {
        local a`i'_`x' = r(`x'_`i')
    }
}
assert "`wv_'" == "_ipw_pswt"
assert `nl_' == 2
assert `lev1' == 0 & `lev2' == 1
di as res "1. run ok; r(weightvar) is the truncated weight"

* ---- 2. truncation bounds, shares, truncated weight ---------------------------------------------------------------------------------
assert reldif(lo_s, 0.0092329060) < 1e-5
assert reldif(hi_s, 6.0235830652) < 1e-5
* definition of the person-time weighted lower quantile: share(w <= q) >= p, share(w < q) < p
gen double len_ = t1 - t0
gen byte dg_ = (adh == 1) & !missing(_ipw_psw) & t1 > t0
quietly summarize len_ if dg_, meanonly
local tot = r(sum)
quietly summarize len_ if dg_ & _ipw_psw <= lo_s, meanonly
assert r(sum) / `tot' >= 0.01 * (1 - 1e-9)
quietly summarize len_ if dg_ & _ipw_psw <  lo_s, meanonly
assert r(sum) / `tot' <  0.01
local sl = r(sum) / `tot'
quietly summarize len_ if dg_ & _ipw_psw <= hi_s, meanonly
assert r(sum) / `tot' >= 0.99 * (1 - 1e-9)
quietly summarize len_ if dg_ & _ipw_psw <  hi_s, meanonly
assert r(sum) / `tot' <  0.99
assert abs(`slo_' - `sl') < 1e-12
assert abs(`slo_' - 0.00997) < 1e-4 & abs(`shi_' - 0.00998) < 1e-4
* truncated weight: equal inside the bounds, set to the bounds outside, missing exactly where psw is missing
assert missing(_ipw_pswt) == missing(_ipw_psw)
assert reldif(_ipw_pswt, min(max(_ipw_psw, lo_s), hi_s)) < 1e-14 if !missing(_ipw_psw)
assert _ipw_pswt >= lo_s & _ipw_pswt <= hi_s if !missing(_ipw_pswt)
assert _ipw_pswt == _ipw_psw if _ipw_psw >= lo_s & _ipw_psw <= hi_s & !missing(_ipw_psw)
di as res "2. person-time weighted percentiles, shares and truncated weight: ok"

* ---- 3. diagnostics (person-time weighted), overall and by arm, against the Python reference ------------------------------------------
assert `all_n_used' == 75182
assert reldif(`all_pt', 118443.1655) < 1e-6 & reldif(`all_mean', 0.97757218) < 1e-5 & reldif(`all_min', 0.00085176) < 1e-5
assert reldif(`all_p1', 0.00923291) < 1e-5 & reldif(`all_p50', 0.81501708) < 1e-5 & reldif(`all_p99', 6.02358307) < 1e-5
assert reldif(`all_max', 15.82371047) < 1e-5 & reldif(`all_ess', 0.45611578) < 1e-5
assert `a1_n' == 37043 & reldif(`a1_pt', 58319.5651) < 1e-6 & reldif(`a1_mean', 1.12640776) < 1e-5
assert reldif(`a1_p1', 0.00634033) < 1e-5 & reldif(`a1_p50', 0.88543061) < 1e-5 & reldif(`a1_p99', 6.54839498) < 1e-5
assert reldif(`a1_max', 14.06583325) < 1e-5 & reldif(`a1_ess', 0.47700588) < 1e-5
assert `a2_n' == 38139 & reldif(`a2_pt', 60123.6005) < 1e-6 & reldif(`a2_mean', 0.83320247) < 1e-5
assert reldif(`a2_p1', 0.02272328) < 1e-5 & reldif(`a2_p50', 0.76110270) < 1e-5 & reldif(`a2_p99', 4.57495231) < 1e-5
assert reldif(`a2_max', 15.82371047) < 1e-5 & reldif(`a2_ess', 0.44863962) < 1e-5
assert abs((`a1_pt' + `a2_pt') - `all_pt') < 1e-6
di as res "3. diagnostics overall and by arm match the Python reference: ok"

* ---- 4. on-treatment Cox with the truncated weight --------------------------------------------------------------------------------------
stset t1 [pw=_ipw_pswt] if adh == 1, time0(t0) failure(inf)
stcox i.arm, breslow vce(cluster id)
di "truncated p1/p99:  coef " %12.8f _b[1.arm] "  SE " %12.8f _se[1.arm]
assert abs(_b[1.arm]  - (-0.4586770716)) < 2e-5
assert abs(_se[1.arm] - 0.1809860006)    < 2e-5
di as res "4. Cox on truncated weights matches the Python reference: ok"

* ---- 5. without truncate()/arm(): no truncated variable, weightvar is the untruncated weight; row order does not matter ---------------------------
drop _ipw_pswt
ipwadh_eval, time0(t0) time1(t1) adherence(adh) replace
local wv2_ = "`r(weightvar)'"
assert "`wv2_'" == "_ipw_psw"
capture confirm variable _ipw_pswt
assert _rc != 0
preserve
set seed 20261010
gen double u_ = runiform()
sort u_
drop u_
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) truncate(1 99) replace
assert abs(r(trunc_lo) - lo_s) < 1e-12 & abs(r(trunc_hi) - hi_s) < 1e-12
restore
di as res "5. no-truncation call and row-order invariance of the bounds: ok"

* ---- 6. input errors ----------------------------------------------------------------------------------------------------------------------
capture ipwadh_eval, time0(t0) time1(t1) truncate(99 1) replace
assert _rc == 198
capture ipwadh_eval, time0(t0) time1(t1) truncate(1 99)                // _ipw_psw exists, no -replace-
assert _rc == 110
di as res "6. input errors (198, 110): ok"

* ---- 7. bands: truncation percentiles are person-time weighted, so they hardly depend on the split (informative, no assertion) ---------------
use "`d'/episodes.dta", clear
ipwadhc, `base' `lists' parts replace
ipwadh_split, `sp' bandpoints(10)
ipwadh_eval, time0(t0) time1(t1) evalat(end) adherence(adh) arm(arm) truncate(1 99) replace
di as txt "bandpoints(10), evalat(end): bounds " %9.6f r(trunc_lo) " / " %9.6f r(trunc_hi) "   (exact split: " %9.6f lo_s " / " %9.6f hi_s ")"

di as res "test_ipwadh_eval_trunc.do: all assertions passed"
