* run_ipwadh_mi_fullsize.do -- FULL-SIZE run: S8_N2000 (69,645 episode rows, ~7.5M pieces after the exact split), M = 10, ipwadh_mi with keep().
* Purpose (open item 2, handoff plan section 7): measure MEMORY and TIME of  ipwadh_mi -> mi stset -> mi estimate: stcox  at full size, and check the pooled result.
* Stata 19.  Run from the folder holding ipwadh.ado, ipwadh_mi.ado (0.2.1), ipwadhc.ado, ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do, ipwadh_eval.ado, and data/S8_N2000/episodes.dta.
* Run it with  do run_ipwadh_mi_fullsize.do  in a FRESH Stata session (NOT through cscript: cscript closes all logs and clears globals).
* Close other big programs first.  Expected: roughly 1 hour in total; RAM: see the notes at the end of this header.
*
* WHAT THIS RUN DECIDES
*   (a) does the FLONG result (M+1 = 11 copies of ~7.5M pieces = ~83M rows, with keep()) fit in this computer's RAM, and how long do ipwadh_mi / mi stset / mi estimate: stcox take?
*   (b) is the pooled result equal to the earlier "within" run (test_ipwadh_S8_mi.do, same data, same imputation seed 20261010, same imputation model)?
*         earlier within-run (Rubin, M=10):  untruncated log HR -0.39316  SE 0.05414      truncated p1/p99 -0.39035  SE 0.05031
*         complete-data benchmark (S3):      untruncated      -0.391823   SE 0.053725      truncated        -0.390331  SE 0.050138
*       Same imputations => the mi estimate numbers should match the within-run numbers up to rounding (stcox i.arm, Breslow, cluster id).
*       If the imputations differ (different Stata update / seed handling) the soft checks against the complete-data benchmark still apply.
*
* IF IT RUNS OUT OF MEMORY (r(901) / r(909) / r(900) / Windows paging to a crawl):  set  local mistyle wide  (cheaper; table only, no mi stset),
*   or keep flong but a smaller keep() list (keep() already drops everything not needed). Report the error and the output of -memory- (printed below).
*
* MEMORY ARITHMETIC (rough): flong keeps ~12 variables + 4 weights, ~140-170 bytes/row x 83M rows = ~12-14 GB, plus mi stset variables (~25 bytes/row) and the temporary
*   frames of ipwadh_mi (the unsplit mi data, the 7.5M-piece work frame ~0.5-1 GB, Mata matrices ~0.8 GB).  If the computer has < 16 GB RAM expect trouble: use  local mistyle wide.

clear all
set more off
capture log close _all
log using "ipwadh_fullsize_run.log", replace text
capture noisily set max_memory .          // a finite max_memory makes Stata itself refuse (lesson of the force tests); this is allowed on Windows

local M 10
local seed 20261010
local mistyle flong          // flong (needed for mi stset / mi estimate: stcox) | wide (weights tables only)
local run_cov 0              // 1 = ALSO fit stcox i.arm c.age i.male (the model that contains exactly the stcoxvars; ~2 x 15 more minutes)
local keepvars "age male"

global nfail = 0
global nchk = 0
capture program drop chk
program define chk
    args label value expected tol mode
    if "`mode'" == "" local mode "abs"
    if "`mode'" == "abs"  local dif = abs(`value' - `expected')
    if "`mode'" == "rel2" local dif = abs(`value' - `expected') / abs(`expected')
    global nchk = $nchk + 1
    if `dif' < `tol' {
        di as txt "  PASS  `label': " %12.6f `value' "  (reference " %12.6f `expected' ")"
    }
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label': " %12.6f `value' "  (reference " %12.6f `expected' ", `mode' diff " %10.3e `dif' ")"
    }
end

di as res "run_ipwadh_mi_fullsize.do   started " c(current_date) " " c(current_time)
di as txt "Stata " c(stata_version) "  " c(os) "  edition " c(edition) "  MP cores " c(processors) "  max_memory " c(max_memory)
query memory
foreach c in ipwadh ipwadh_mi ipwadhc ipwadh_split ipwadh_eval ipwadh_covlists {
    capture noisily which `c'
}

local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local opts "`base' failure(inf) `lists'"

timer clear

* ===== 0. data + imputation (EXACTLY as test_ipwadh_S8_mi.do, so the imputations are the same) =============================================================
timer on 9
capture confirm file "S8_N2000_mi_M`M'.dta"
if _rc {
    use "data/S8_N2000/episodes.dta", clear
    assert _N == 69645
    stset t1, time0(t0) failure(inf)
    quietly sts generate na = na
    stset, clear
    quietly replace na = 0 if missing(na)
    mi set wide
    mi register imputed biomarker abx
    mi register regular arm age male adh tsi ninf inf dead ltfu na t0 t1
    mi impute chained (regress) biomarker (logit, augment) abx = arm age male adh i.tsi ninf inf dead ltfu na t1, add(`M') rseed(`seed')
    save "S8_N2000_mi_M`M'.dta", replace
    di as res "0. imputation done and saved: S8_N2000_mi_M`M'.dta"
}
else di as res "0. using the saved imputed data S8_N2000_mi_M`M'.dta (delete the file to re-impute)"
use "S8_N2000_mi_M`M'.dta", clear
quietly mi describe
assert r(M) == `M'
timer off 9

* ===== A. ipwadh_mi: weights in every imputation, ONE exact split, result frame =============================================================================
di as res _n "=== A. ipwadh_mi (mistyle(`mistyle'), M = `M', keep(`keepvars')) ===   " c(current_time)
timer on 1
capture noisily ipwadh_mi, `opts' frame(res, replace) keep(`keepvars') mistyle(`mistyle')
local rcA = _rc
timer off 1
if `rcA' {
    di as err _n "ipwadh_mi FAILED with rc = `rcA'.  Send me this log (header, -query memory- and the error above)."
    di as err "If rc is 901/909/900 (memory): set  local mistyle wide  and re-run; the imputed data are saved, so step 0 is skipped."
    timer list
    log close
    exit `rcA'
}
frame change res
di as txt _n "result frame after ipwadh_mi:"
describe, short
memory
quietly count if _mi_m == 0
local Np = r(N)
di as txt "pieces per imputation: " %12.0fc `Np'

if "`mistyle'" != "flong" {
    di as res _n "mistyle(wide): only the weights table is available (no mi stset). Weighted baseline table, adherent pieces:"
    timer on 2
    mi estimate: mean age biomarker [pw=_ipw_ptw] if adh == 1, over(arm)
    timer off 2
    timer list
    di as txt "Done (wide). Send me the whole log."
    log close
    exit
}

* ===== B. mi stset + mi estimate: stcox, untruncated weights ================================================================================================
di as res _n "=== B. mi stset [pw=_ipw_psw] + mi estimate: stcox i.arm, breslow vce(cluster id) ===   " c(current_time)
timer on 2
mi stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
timer off 2
memory
timer on 3
mi estimate, post: stcox i.arm, breslow vce(cluster id)
timer off 3
scalar ru_q = _b[1.arm]
scalar ru_se = _se[1.arm]
di as txt "  untruncated: log HR " %9.5f scalar(ru_q) "   SE " %8.5f scalar(ru_se) "   HR " %7.4f exp(scalar(ru_q))

* ===== C. truncated weights (p1/p99) ================================================================================================================================
di as res _n "=== C. mi stset [pw=_ipw_pswt] + mi estimate: stcox i.arm ===   " c(current_time)
timer on 4
capture noisily mi stset t1 [pw=_ipw_pswt] if adh == 1, time0(t0) failure(inf)
if _rc {
    di as txt "  (repeated mi stset refused; trying mi stset, clear first)"
    capture noisily mi stset, clear
    mi stset t1 [pw=_ipw_pswt] if adh == 1, time0(t0) failure(inf)
}
timer off 4
timer on 5
mi estimate, post: stcox i.arm, breslow vce(cluster id)
timer off 5
scalar rt_q = _b[1.arm]
scalar rt_se = _se[1.arm]
di as txt "  truncated p1/p99: log HR " %9.5f scalar(rt_q) "   SE " %8.5f scalar(rt_se) "   HR " %7.4f exp(scalar(rt_q))

* ===== D. optional: the model that contains exactly the stcoxvars =================================================================================================
if `run_cov' == 1 {
    di as res _n "=== D. mi estimate: stcox i.arm c.age i.male (exactly the stcoxvars; truncated, then untruncated) ===   " c(current_time)
    timer on 6
    mi estimate, post: stcox i.arm c.age i.male, breslow vce(cluster id)
    scalar dt_q = _b[1.arm]
    scalar dt_se = _se[1.arm]
    capture noisily mi stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
    if _rc {
        capture noisily mi stset, clear
        mi stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
    }
    mi estimate, post: stcox i.arm c.age i.male, breslow vce(cluster id)
    scalar du_q = _b[1.arm]
    scalar du_se = _se[1.arm]
    timer off 6
    di as txt "  with age and male: untruncated log HR " %9.5f scalar(du_q) " SE " %8.5f scalar(du_se) "   truncated log HR " %9.5f scalar(dt_q) " SE " %8.5f scalar(dt_se)
}

* ===== E. checks ======================================================================================================================================================
di as res _n "=== E. checks ==="
di as txt "-- against the earlier within-imputation run (same data, same imputations expected; rounded to 5 decimals there) --"
chk "log HR untruncated vs within run"  scalar(ru_q)  -0.39316 5e-5
chk "SE untruncated vs within run"      scalar(ru_se)  0.05414 5e-5
chk "log HR truncated vs within run"    scalar(rt_q)  -0.39035 5e-5
chk "SE truncated vs within run"        scalar(rt_se)  0.05031 5e-5
di as txt "-- against the complete-data benchmark (soft: tolerance as in test_ipwadh_S8_mi.do) --"
chk "log HR untruncated vs complete data" scalar(ru_q) -0.391823 0.04
chk "SE untruncated vs complete data"     scalar(ru_se) 0.053725 0.20 rel2
chk "log HR truncated vs complete data"   scalar(rt_q) -0.390331 0.04
chk "SE truncated vs complete data"       scalar(rt_se) 0.050138 0.20 rel2
di as txt "(a FAIL in the first four checks with passes in the last four means: different imputations, not a defect; send me the log)"

* ===== F. timing summary ==============================================================================================================================================
di as res _n "=== F. time and memory ==="
timer list
di as txt "timer 9 = imputation (0 if reused), 1 = ipwadh_mi, 2 = mi stset (untrunc), 3 = mi estimate: stcox untrunc, 4 = mi stset (trunc), 5 = mi estimate: stcox trunc, 6 = optional model with age male"
di as txt "(seconds).  Reference at N=2000, ONE imputation: ipwadhc 4.8 s, split 12.5 s, eval 16.4 s, stset + ONE stcox 71.8 s."
memory
di as res _n "run_ipwadh_mi_fullsize.do: $nchk checks, $nfail failed.   finished " c(current_date) " " c(current_time)
di as txt "Send me the whole log (ipwadh_fullsize_run.log)."
frame change default
log close
