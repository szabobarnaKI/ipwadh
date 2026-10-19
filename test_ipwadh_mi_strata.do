* test_ipwadh_mi_strata.do -- certification of  mi estimate: stcox ..., strata()  on the result of ipwadh_mi / ipwadh (flong, v0.2.1).  NOT YET RUN IN STATA.  Stata 16+ (run in Stata 19).
* Run from the folder containing ipwadh.ado (0.2.1), ipwadh_mi.ado (0.2.1), ipwadhc.ado, ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do, ipwadh_eval.ado and data/S8_N2000/episodes.dta.
* Data: the first 300 persons of S8_N2000, M = 3 imputations (as test_ipwadh_mi_frame.do).  Part of the FULL set of cscript_ipwadh.do.
* RULES under test (from the design of the stratified split, LOG.md / handoff plan section 3.2)
*   (i)  ipwadh_mi/ipwadh strata(varlist) cuts the rows only at the failure times of the row's own stratum.  The stratified split is exact ONLY for a Cox model stratified by the same variables.
*   (ii) the Cox model must be arm + exactly the stcoxvars; a variable used in strata() of stcox is then NOT in stcoxvars (it is a confounder: put it in adhpredictors()).
*   (iii) the strata variables of the split may not be imputed / passive (structure of the data): refused with rc 198.  Stratifying the COX model by an imputed variable is fine when the split is not stratified.
* What it certifies (every comparison is against an INDEPENDENT reference loop:  mi extract m -> ipwadhc parts -> ipwadh_split [strata()] -> ipwadh_eval -> stset [pw] -> stcox ..., strata(); pooled by Rubin's rules in Mata)
*   1  error paths: strata() with an imputed variable -> 198 (ipwadh_mi and ipwadh), no frame left behind
*   2  A  split strata(male),       Cox strata(male)         : pieces, weights, pooled b, total variance (Rubin); also with the TRUNCATED weight (mistset pswt)
*   3  B  split not stratified,     Cox strata(male)         : same checks; pooled b equals case A (exactness of the stratified split)
*   4  C  split strata(tsi) (time-varying stratum, constant within a row),  Cox strata(tsi)
*   5  D  split strata(male tsi) (two variables),  Cox strata(male tsi)
*   6  E  split not stratified,     Cox strata(abx) with abx IMPUTED (the stratum differs between imputations)
*   7  wrapper:  ipwadh ... strata(male) mistset(psw) == manual mi stset on ipwadh_mi; the stratified split has fewer pieces than the exact unstratified one
*   8  INFORMATION only (not counted): UNSTRATIFIED Cox on the stratified split is not exact -> the difference to the exact result is printed
clear all
set more off
global nfail = 0
global nchk = 0
local M 3
local seed 20261010

capture program drop chk
program define chk
    args label value expected tol
    local dif = abs(`value' - `expected')
    global nchk = $nchk + 1
    if `dif' < `tol' {
        di as txt "  PASS  `label': " %14.8g `value' "  (expected " %14.8g `expected' ")"
    }
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label': " %14.8g `value' "  (expected " %14.8g `expected' ", diff " %10.3e `dif' ")"
    }
end

capture program drop chkrc
program define chkrc
    args label rc expected
    global nchk = $nchk + 1
    if `rc' == `expected' di as txt "  PASS  `label': rc = `rc'"
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label': rc = `rc', expected `expected'"
    }
end

capture program drop chktrue
program define chktrue
    args label cond
    global nchk = $nchk + 1
    if `cond' di as txt "  PASS  `label'"
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label'"
    }
end

* data and imputation (as in test_ipwadh_mi_frame.do)
capture program drop mkmi
program define mkmi
    args style M seed
    use "data/S8_N2000/episodes.dta", clear
    keep if id <= 300
    quietly stset t1, time0(t0) failure(inf)
    quietly sts generate na = na
    stset, clear
    quietly replace na = 0 if missing(na)
    mi set wide
    mi register imputed biomarker abx
    mi register regular arm age male adh tsi ninf inf dead ltfu na t0 t1
    quietly mi impute chained (regress) biomarker (logit, augment) abx = arm age male adh i.tsi ninf inf dead ltfu na t1, add(`M') rseed(`seed')
    if "`style'" != "wide" quietly mi convert `style', clear
end

* compare the weights / imputed covariates of a flong result frame with the reference matrices REFK, REFW1..REFWM (from refloop)
capture program drop cmpw
program define cmpw
    args fr label M Np
    frame `fr' {
        quietly count
        chk "`label': rows of the result" r(N) `=`Np' * (`M' + 1)' 0.5
        mata: RA = st_data(., ("_ipw_psw", "_ipw_pswt", "_ipw_ptw", "_ipw_ptwt", "biomarker", "abx"))
        mata: KA = st_data(., "id t0 t1")
        forvalues m = 1/`M' {
            mata: RW = RA[((`m' * `Np' + 1)::((`m' + 1) * `Np')), .]
            mata: st_numscalar("dkey", max(abs(KA[((`m' * `Np' + 1)::((`m' + 1) * `Np')), .] - REFK)))
            mata: st_numscalar("dmiss", sum((RW :>= .) :!= (REFW`m' :>= .)))
            mata: st_numscalar("dval", max(abs(editmissing(RW, 0) - editmissing(REFW`m', 0))))
            chk "`label' m=`m': piece keys (id t0 t1) equal the reference" scalar(dkey) 0 1e-12
            chk "`label' m=`m': cells whose missing status differs" scalar(dmiss) 0 0.5
            chk "`label' m=`m': max abs diff psw/pswt/ptw/ptwt/biomarker/abx" scalar(dval) 0 1e-9
        }
        capture mata: mata drop RA KA RW
    }
end

* reference loop.  Inputs (globals): LISTS (weight-model lists), SPLITST (strata() of the split, may be empty), COXSPEC (stcox command after "stcox "), WVAR (weight variable of the Cox model).
* Outputs: Mata REFK, REFW1..REFWM; matrices BM (M x k per-imputation coefficients) and WM (mean of the per-imputation variances); global NREF (pieces)
capture program drop refloop
program define refloop
    args M
    local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu)"
    local sp ""
    if "$SPLITST" != "" local sp "strata($SPLITST)"
    capture matrix drop BM
    capture matrix drop WM
    forvalues m = 1/`M' {
        capture frame drop Fr
        frame copy default Fr
        frame Fr {
            quietly mi extract `m', clear
            quietly ipwadhc, `base' $LISTS parts replace
            quietly ipwadh_split, id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu) `sp'
            quietly ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) truncate(1 99) replace
            quietly gen double _ipw_ptw = _ipw_psw * (t1 - t0)
            quietly gen double _ipw_ptwt = _ipw_pswt * (t1 - t0)
            mata: REFK = st_data(., "id t0 t1")
            mata: REFW`m' = st_data(., "_ipw_psw _ipw_pswt _ipw_ptw _ipw_ptwt biomarker abx")
            quietly stset t1 [pw=$WVAR] if adh == 1, time0(t0) failure(inf)
            quietly stcox $COXSPEC
            matrix bb = e(b)
            matrix vv = e(V)
            local nref = _N
        }
        frame drop Fr
        matrix BM = (nullmat(BM) \ bb)
        if `m' == 1 matrix WM = vv / `M'
        else matrix WM = WM + vv / `M'
    }
    global NREF = `nref'
end

* mi stset + mi estimate: stcox on the result frame $RES with $COXSPEC and the weight $WVAR; checks against the reference (BM, WM); leaves matrices bmi, vmi
capture program drop docase
program define docase
    args label tol
    capture matrix drop bmi
    capture matrix drop vmi
    capture matrix drop wmi
    frame $RES {
        capture noisily mi stset t1 [pw=$WVAR] if adh == 1, time0(t0) failure(inf)
        local rc_s = _rc
        capture noisily mi estimate, post: stcox $COXSPEC
        local rc_h = _rc
        if `rc_h' == 0 {
            matrix bmi = e(b)
            matrix vmi = e(V)
            capture matrix wmi = e(W_mi)
        }
    }
    chkrc "`label': mi stset [pw=$WVAR]" `rc_s' 0
    chkrc "`label': mi estimate: stcox $COXSPEC" `rc_h' 0
    if `rc_h' == 0 {
        mata: BMr = st_matrix("BM")
        mata: bbar = mean(BMr)
        mata: Bb = variance(BMr)
        mata: Tt = st_matrix("WM") + (1 + 1 / rows(BMr)) * Bb
        mata: st_numscalar("d_b", max(abs(st_matrix("bmi") - bbar)))
        mata: st_numscalar("d_v", max(abs(st_matrix("vmi") - Tt)))
        chk "`label': pooled coefficients = mean over m of per-m stset/stcox" scalar(d_b) 0 `tol'
        chk "`label': total variance = W + (1+1/M) B from the per-m fits" scalar(d_v) 0 `tol'
        capture confirm matrix wmi
        if _rc == 0 {
            mata: st_numscalar("d_w", max(abs(st_matrix("wmi") - st_matrix("WM"))))
            chk "`label': within-imputation variance = mean of per-m V" scalar(d_w) 0 `tol'
        }
        else di as txt "  (e(W_mi) not available: within-imputation variance compared through the total variance only)"
        mata: st_numscalar("nzero", sum(st_matrix("bmi") :>= .))
        chk "`label': no missing pooled coefficient" scalar(nzero) 0 0.5
        di as txt "  pooled coefficients:"
        matrix list bmi, format(%9.5f)
        capture mata: mata drop BMr bbar Bb Tt
    }
end

local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu) failure(inf)"
global LISTS1 "stcoxvars(c.age) adhpredictors(i.male i.tsi i.abx c.biomarker)"
global LISTS2 "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local tolcase 1e-8

mkmi wide `M' `seed'
quietly datasignature
local sig0 "`r(datasignature)'"

* ===== 1. error paths ==============================================================================================================================
di as res _n "=== 1. error paths ==="
capture noisily ipwadh_mi, `base' $LISTS1 strata(abx) frame(s_err)
chkrc "ipwadh_mi: strata() with an imputed variable" `=_rc' 198
capture noisily ipwadh, `base' $LISTS1 strata(abx) frame(s_err)
chkrc "ipwadh: strata() with an imputed variable" `=_rc' 198
capture frame s_err: describe
global nchk = $nchk + 1
if _rc == 0 {
    global nfail = $nfail + 1
    di as err "  FAIL  a result frame s_err was left behind by a failed call"
}
else di as txt "  PASS  no frame left behind by failed calls"

* ===== 2. case A: split strata(male), Cox strata(male) =============================================================================================
di as res _n "=== 2. A: split strata(male), mi estimate: stcox, strata(male) ==="
ipwadh_mi, `base' $LISTS1 strata(male) frame(rA)
frame rA: quietly count if _mi_m == 0
local npA = r(N)
di as txt "  pieces per imputation (A): `npA'"
global LISTS "$LISTS1"
global SPLITST "male"
global COXSPEC "i.arm c.age, strata(male) breslow vce(cluster id)"
global WVAR "_ipw_psw"
global RES "rA"
refloop `M'
chk "A: pieces of the result = reference" `npA' $NREF 0.5
cmpw rA "A" `M' `npA'
docase "A (untruncated)" `tolcase'
matrix bA = bmi

global WVAR "_ipw_pswt"
refloop `M'
docase "A (truncated pswt)" `tolcase'

* ===== 3. case B: split NOT stratified, Cox strata(male) ===========================================================================================
di as res _n "=== 3. B: exact unstratified split, mi estimate: stcox, strata(male) ==="
ipwadh_mi, `base' $LISTS1 frame(rB)
frame rB: quietly count if _mi_m == 0
local npB = r(N)
di as txt "  pieces per imputation (B): `npB'  (A: `npA')"
chktrue "the stratified split has fewer pieces than the exact unstratified split" "`npA' < `npB'"
global SPLITST ""
global WVAR "_ipw_psw"
global RES "rB"
refloop `M'
chk "B: pieces of the result = reference" `npB' $NREF 0.5
cmpw rB "B" `M' `npB'
docase "B" `tolcase'
matrix bB = bmi
mata: st_numscalar("dAB", max(abs(st_matrix("bA") - st_matrix("bB"))))
chk "A = B: stratified split gives the same stratified Cox estimate as the exact split" scalar(dAB) 0 1e-5

* ===== 4. case C: time-varying stratum tsi ===========================================================================================================
di as res _n "=== 4. C: split strata(tsi), mi estimate: stcox, strata(tsi) ==="
ipwadh_mi, `base' $LISTS2 strata(tsi) frame(rC)
frame rC: quietly count if _mi_m == 0
local npC = r(N)
di as txt "  pieces per imputation (C): `npC'"
global LISTS "$LISTS2"
global SPLITST "tsi"
global COXSPEC "i.arm c.age i.male, strata(tsi) breslow vce(cluster id)"
global RES "rC"
refloop `M'
chk "C: pieces of the result = reference" `npC' $NREF 0.5
cmpw rC "C" `M' `npC'
docase "C" `tolcase'

* ===== 5. case D: two strata variables ===============================================================================================================
di as res _n "=== 5. D: split strata(male tsi), mi estimate: stcox, strata(male tsi) ==="
ipwadh_mi, `base' $LISTS1 strata(male tsi) frame(rD)
frame rD: quietly count if _mi_m == 0
local npD = r(N)
di as txt "  pieces per imputation (D): `npD'  (A: `npA')"
chktrue "two strata variables: not more pieces than strata(male)" "`npD' <= `npA'"
global LISTS "$LISTS1"
global SPLITST "male tsi"
global COXSPEC "i.arm c.age, strata(male tsi) breslow vce(cluster id)"
global RES "rD"
refloop `M'
chk "D: pieces of the result = reference" `npD' $NREF 0.5
cmpw rD "D" `M' `npD'
docase "D" `tolcase'

* ===== 6. case E: imputed stratum in the Cox model (split not stratified) ============================================================================
di as res _n "=== 6. E: exact split, mi estimate: stcox, strata(abx) with abx imputed ==="
ipwadh_mi, `base' $LISTS2 frame(rE)
frame rE: quietly count if _mi_m == 0
local npE = r(N)
global LISTS "$LISTS2"
global SPLITST ""
global COXSPEC "i.arm c.age i.male, strata(abx) breslow vce(cluster id)"
global RES "rE"
refloop `M'
chk "E: pieces of the result = reference" `npE' $NREF 0.5
cmpw rE "E" `M' `npE'
docase "E" `tolcase'

* ===== 7. wrapper =====================================================================================================================================
di as res _n "=== 7. ipwadh ... strata(male) mistset(psw) ==="
ipwadh, `base' $LISTS1 strata(male) frame(rW) mistset(psw)
frame rW {
    quietly count if _mi_m == 0
    local npW = r(N)
    capture noisily mi estimate, post: stcox i.arm c.age, strata(male) breslow vce(cluster id)
    local rcW = _rc
    matrix bW = e(b)
}
chk "wrapper: pieces = ipwadh_mi" `npW' `npA' 0.5
chkrc "wrapper: mi estimate: stcox, strata(male) after mistset(psw)" `rcW' 0
mata: st_numscalar("dW", max(abs(st_matrix("bW") - st_matrix("bA"))))
chk "wrapper + mistset(psw) = ipwadh_mi + manual mi stset (case A)" scalar(dW) 0 1e-10

* ===== 8. information only: UNSTRATIFIED Cox on the stratified split ================================================================================
di as res _n "=== 8. information (not counted): unstratified stcox on the strata(male) split vs the exact split ==="
foreach f in rA rB {
    frame `f' {
        quietly mi stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
        quietly mi estimate, post: stcox i.arm c.age, breslow vce(cluster id)
        matrix bU`f' = e(b)
    }
}
mata: st_numscalar("dU", max(abs(st_matrix("bUrA") - st_matrix("bUrB"))))
di as txt "  max abs difference of the pooled coefficients, unstratified Cox on split A vs on the exact split B: " as res %10.3e scalar(dU)
di as txt "  (an unstratified Cox model is only exact on the exact UNSTRATIFIED split; use strata() in the split only together with the same strata in stcox)"

foreach f in rA rB rC rD rE rW {
    capture frame drop `f'
}
quietly datasignature
global nchk = $nchk + 1
if "`r(datasignature)'" == "`sig0'" di as txt "  PASS  datasignature of the user's mi data unchanged"
else {
    global nfail = $nfail + 1
    di as err "  FAIL  the user's mi data changed"
}

di as res _n "test_ipwadh_mi_strata.do: $nchk checks, $nfail failed"
di as txt "Send me the whole log."
