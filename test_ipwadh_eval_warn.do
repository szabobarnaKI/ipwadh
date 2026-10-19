* test_ipwadh_eval_warn.do -- run in Stata 19 from the folder containing ipwadhc.ado (v0.3.0), ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do,
* ipwadh_eval.ado (v0.3.0) and data/.   NOT YET RUN.   Data: S3_N200 (seed 1), exact split (77,718 pieces), evalat(end), death model only (as test_ipwadh_eval_trunc.do).
* Tests the weight warning of ipwadh_eval 0.3.0: esscaution() esssevere() maxwarn() nowarn, returns r(ess_min) r(warn) r(warn_max) r(ess_trunc).
* Expected numbers (independent Python reference, ess_ref.py / ref_core.py):  ESS overall 0.45611578, arm0 0.47700588, arm1 0.44863962, max weight 15.82371047,
*   ESS of the weight truncated at person-time p1/p99 = 0.5100566993.
* ALSO re-run test_ipwadh_eval.do and test_ipwadh_eval_trunc.do (ipwadh_eval 0.3.0 must not change any of their results).
clear all
set more off
global nfail = 0
global nchk = 0
local d "data/S3_N200"
local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local sp "id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu)"

capture program drop chk
program define chk
    args label value expected tol
    global nchk = $nchk + 1
    if abs(`value' - `expected') < `tol' di as txt "  PASS  `label': " %14.8f `value' "  (expected " %14.8f `expected' ")"
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label': " %14.8f `value' "  (expected " %14.8f `expected' ")"
    }
end

use "`d'/episodes.dta", clear
ipwadhc, `base' `lists' parts replace
ipwadh_split, `sp'
assert _N == 77718

* ===== 1. defaults: ESS 45-50% is stable -> no warning; r(ess_min) is the smallest of overall and by-arm ESS ========================================
di as res _n "=== 1. defaults, with arm() and truncate(1 99) ==="
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) truncate(1 99) replace
chk "ess_min (= arm 1)" r(ess_min) 0.44863962 1e-6
chk "ess (overall)" r(ess) 0.45611578 1e-6
chk "warn (none expected)" r(warn) 0 0.5
chk "warn_max (max 15.8 < 100)" r(warn_max) 0 0.5
chk "ess_trunc" r(ess_trunc) 0.5100566993 1e-6
chk "trunc_hi unchanged" r(trunc_hi) 6.0235830652 1e-6

* ===== 2. thresholds ====================================================================================================================================
di as res _n "=== 2. esscaution() / esssevere() ==="
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) esscaution(50) replace
chk "esscaution(50): caution" r(warn) 1 0.5
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) esscaution(60) esssevere(50) replace
chk "esscaution(60) esssevere(50): severe" r(warn) 2 0.5
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) esscaution(30) esssevere(20) replace
chk "esscaution(30) esssevere(20): none" r(warn) 0 0.5
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) esscaution(45) esssevere(44) replace
chk "ess_min 44.86 between 44 and 45: caution" r(warn) 1 0.5
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) esscaution(100) esssevere(44.9) replace
chk "ess_min 44.86 < 44.9: severe" r(warn) 2 0.5

* ===== 3. the by-arm minimum is used when arm() is given, the overall ESS otherwise ============================================================
di as res _n "=== 3. arm() versus overall ==="
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) esscaution(45) esssevere(10) replace
chk "with arm(): ess_min 44.86 < 45 -> caution" r(warn) 1 0.5
ipwadh_eval, time0(t0) time1(t1) adherence(adh) esscaution(45) esssevere(10) replace
chk "without arm(): ess_min = overall 45.61 > 45 -> none" r(warn) 0 0.5
chk "without arm(): r(ess_min)" r(ess_min) 0.45611578 1e-6
assert missing(r(ess_trunc))
di as txt "  PASS  no truncate(): r(ess_trunc) is missing"

* ===== 4. maxwarn() =========================================================================================================================================
di as res _n "=== 4. maxwarn() ==="
ipwadh_eval, time0(t0) time1(t1) adherence(adh) maxwarn(10) replace
chk "maxwarn(10): max 15.8 > 10" r(warn_max) 1 0.5
ipwadh_eval, time0(t0) time1(t1) adherence(adh) maxwarn(15.9) replace
chk "maxwarn(15.9): max 15.82 < 15.9" r(warn_max) 0 0.5

* ===== 5. nowarn: returns are unchanged (the printed text is not tested, look at the output of section 2 for the messages) ==========================
di as res _n "=== 5. nowarn ==="
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) esscaution(60) esssevere(50) maxwarn(10) nowarn replace
chk "nowarn: warn still returned" r(warn) 2 0.5
chk "nowarn: warn_max still returned" r(warn_max) 1 0.5

* ===== 6. input errors ======================================================================================================================================
di as res _n "=== 6. errors ==="
foreach opt in "esssevere(30) esscaution(20)" "esssevere(0)" "esscaution(101)" "maxwarn(0)" "esssevere(20) esscaution(20)" {
    capture ipwadh_eval, time0(t0) time1(t1) adherence(adh) `opt' replace
    chk "`opt': rc 198" _rc 198 0.5
}

* ===== 7. optional, informational: S6 (near-violation of positivity, session 4g: ESS about 5.7%, max weight about 220) ===============================
local runS6 0
if `runS6' {
    use "data/S6_N2000/episodes.dta", clear
    ipwadhc, `base' `lists' parts replace
    ipwadh_split, `sp'
    ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) truncate(1 99) replace
    di as txt "S6: ess_min " %7.4f r(ess_min) "  warn " r(warn) "  warn_max " r(warn_max) "  ess_trunc " %7.4f r(ess_trunc)
}

di as res _n "test_ipwadh_eval_warn.do: $nchk checks, $nfail failed"
