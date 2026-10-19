* test_ipwadh_mi_frame.do -- certification of ipwadh_mi.ado v0.2.1 (run 2026-10-11 02:14: 99 of 100 PASS; the one FAIL was a wrong expectation of the test, fixed here).  Stata 16+ (run in Stata 19).
* Run from the folder containing ipwadh_mi.ado, ipwadhc.ado (v0.3.0), ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do, ipwadh_eval.ado (v0.2.0)
* and data/S8_N2000/episodes.dta.  Data: the first 300 persons of S8_N2000 (id <= 300), M = 3 imputations.
* v0.2 default result = FLONG mi data: weights are UNREGISTERED super-varying variables (m=0 value 1), so that
*     mi stset t1 [pw=_ipw_psw] if adh==1, ...    mi estimate: stcox ...    mi estimate: mean x [pw=_ipw_ptw], over(arm)      all work natively.
* What it certifies
*   1  error paths (rc 119, 198, 198, 110, no frame left behind), the user's data untouched (datasignature)
*   2  the flong result: row layout (_mi_m, _mi_id), m=0 weights = 1 (ptw = t1-t0), imputed covariates registered imputed, the weights NOT registered
*   3  the weights and imputed covariates of every imputation equal an INDEPENDENT reference loop (mi extract m -> ipwadhc parts -> ipwadh_split -> ipwadh_eval)
*   4  mi estimate: mean [pw=ptw / ptwt] if adh==1, over(arm)  ==  mean over m of per-m weighted means
*   5  mi stset [pw=_ipw_psw] + mi estimate: stcox  ==  mean over m of the per-m stset/stcox (log HRs; also the pooled variance is shown)
*   6  options: frame(, replace), change, notruncate, mistset(psw), keep()
*   7  mistyle(wide) (weights registered imputed; mean only), and flong / mlong SOURCE data
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

capture program drop chktxt
program define chktxt
    args label value expected
    global nchk = $nchk + 1
    if `"`value'"' == `"`expected'"' di as txt "  PASS  `label': `value'"
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label': `value'  (expected `expected')"
    }
end

local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local opts "`base' failure(inf) `lists'"

* ===== 0. data and imputation (as in test_ipwadh_S8_mi.do, 300 persons) ============================================================================
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

* compare the weights / imputed covariates of a result frame with the reference matrices REFW1..REFWM (Mata) and REFK
capture program drop cmpw
program define cmpw
    args fr mstyle label M Np
    frame `fr' {
        quietly count
        local want = cond("`mstyle'" == "flong", `Np' * (`M' + 1), `Np')
        chk "`label': rows of the result" r(N) `want' 0.5
        if "`mstyle'" == "flong" {
            mata: RA = st_data(., ("_ipw_psw", "_ipw_pswt", "_ipw_ptw", "_ipw_ptwt", "biomarker", "abx"))
            mata: KA = st_data(., "id t0 t1")
            mata: st_numscalar("dm0", max(abs(RA[1::`Np', 1..2] :- 1)))
            mata: st_numscalar("dm0b", max(abs(RA[1::`Np', 3..4] - (KA[1::`Np', 3] - KA[1::`Np', 2]) * J(1, 2, 1))))
            chk "`label': m=0 weights psw, pswt = 1" scalar(dm0) 0 1e-12
            chk "`label': m=0 person-time weights = t1 - t0" scalar(dm0b) 0 1e-12
        }
        forvalues m = 1/`M' {
            if "`mstyle'" == "flong" {
                mata: RW = RA[((`m' * `Np' + 1)::((`m' + 1) * `Np')), .]
                mata: st_numscalar("dkey", max(abs(KA[((`m' * `Np' + 1)::((`m' + 1) * `Np')), .] - REFK)))
            }
            else {
                mata: RW = st_data(., ("_`m'__ipw_psw", "_`m'__ipw_pswt", "_`m'__ipw_ptw", "_`m'__ipw_ptwt", "_`m'_biomarker", "_`m'_abx"))
                mata: st_numscalar("dkey", max(abs(st_data(., "id t0 t1") - REFK)))
            }
            mata: st_numscalar("dmiss", sum((RW :>= .) :!= (REFW`m' :>= .)))
            mata: st_numscalar("dval", max(abs(editmissing(RW, 0) - editmissing(REFW`m', 0))))
            chk "`label' m=`m': piece keys (id t0 t1) equal the reference" scalar(dkey) 0 1e-12
            chk "`label' m=`m': cells whose missing status differs" scalar(dmiss) 0 0.5
            chk "`label' m=`m': max abs diff psw/pswt/ptw/ptwt/biomarker/abx" scalar(dval) 0 1e-9
        }
    }
end

mkmi wide `M' `seed'
quietly mi describe
chk "M of the mi data" r(M) `M' 0.5
quietly datasignature
local sig0 "`r(datasignature)'"
local N0 = _N
local style0 "`_dta[_mi_style]'"

* ===== 1. error paths ==============================================================================================================================
di as res _n "=== 1. error paths ==="
capture frame drop plain
frame create plain
frame plain: quietly use "data/S8_N2000/episodes.dta", clear
capture noisily frame plain: ipwadh_mi, `opts' frame(r_plain)
chkrc "non-mi data" `=_rc' 119
capture noisily ipwadh_mi if biomarker > 0, `opts' frame(r_err)
chkrc "if on an imputed variable" `=_rc' 198
capture noisily ipwadh_mi, id(id) time0(t0) time1(t1) adherence(adh) arm(biomarker) failure(inf) dead(dead) ltfu(ltfu) frame(r_err)
chkrc "imputed variable as arm()" `=_rc' 198
capture noisily ipwadh_mi, `opts' frame(r_err) mistyle(long)
chkrc "bad mistyle()" `=_rc' 198
capture noisily ipwadh_mi, `opts' frame(r_err) mistyle(wide) mistset(psw)
chkrc "mistset() with mistyle(wide)" `=_rc' 198
capture frame create r_exist
capture noisily ipwadh_mi, `opts' frame(r_exist)
chkrc "existing frame without replace" `=_rc' 110
capture frame drop r_exist
capture frame r_err: describe
global nchk = $nchk + 1
if _rc == 0 {
    global nfail = $nfail + 1
    di as err "  FAIL  a result frame r_err was left behind by a failed call"
}
else di as txt "  PASS  no frame left behind by failed calls"
capture frame drop plain

* ===== 2. the first successful call (default: flong) ===================================================================================================
di as res _n "=== 2. ipwadh_mi (wide source, default mistyle(flong)) ==="
ipwadh_mi, `opts' frame(res)
chktxt "current frame after the call" "`c(frame)'" "default"
quietly datasignature
chktxt "datasignature of the user's data unchanged" "`r(datasignature)'" "`sig0'"
chk "N of the user's data unchanged" _N `N0' 0.5
chktxt "mi style of the user's data" "`_dta[_mi_style]'" "`style0'"

frame res {
    chktxt "result frame: mi style" "`_dta[_mi_style]'" "flong"
    local Mres = real("`_dta[_mi_M]'")
    chk "result frame: M" `Mres' `M' 0.5
    local iv "`_dta[_mi_ivars]'"
    foreach v in biomarker abx {
        global nchk = $nchk + 1
        if strpos(" `iv' ", " `v' ") == 0 {
            global nfail = $nfail + 1
            di as err "  FAIL  `v' is not registered imputed (ivars: `iv')"
        }
        else di as txt "  PASS  `v' registered imputed"
    }
    local reg "`_dta[_mi_ivars]' `_dta[_mi_pvars]' `_dta[_mi_rvars]'"
    foreach v in _ipw_psw _ipw_pswt _ipw_ptw _ipw_ptwt {
        global nchk = $nchk + 1
        if strpos(" `reg' ", " `v' ") != 0 {
            global nfail = $nfail + 1
            di as err "  FAIL  `v' is registered (it must be unregistered so that it may vary by imputation)"
        }
        else di as txt "  PASS  `v' is unregistered"
    }
    quietly count if _mi_m == 0
    local Np = r(N)
    di as txt "  pieces per imputation: `Np'"
    capture noisily mi varying
}

* ===== 3. the reference loop (independent of ipwadh_mi) =============================================================================================
di as res _n "=== 3. weights per imputation against the reference loop ==="
matrix bref_a = .
matrix bref_t = .
matrix bref_h = .
matrix vref_h = .
forvalues m = 1/`M' {
    capture frame drop Fr
    frame copy default Fr
    frame Fr {
        quietly mi extract `m', clear
        quietly ipwadhc, `base' `lists' parts replace
        quietly ipwadh_split, id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu)
        quietly ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) truncate(1 99) replace
        quietly gen double _ipw_ptw = _ipw_psw * (t1 - t0)
        quietly gen double _ipw_ptwt = _ipw_pswt * (t1 - t0)
        mata: REFK = st_data(., "id t0 t1")
        mata: REFW`m' = st_data(., "_ipw_psw _ipw_pswt _ipw_ptw _ipw_ptwt biomarker abx")
        quietly mean age biomarker [pw=_ipw_ptw] if adh == 1, over(arm)
        matrix ba = e(b)
        quietly mean age biomarker [pw=_ipw_ptwt] if adh == 1, over(arm)
        matrix bt = e(b)
        quietly stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
        quietly stcox i.arm c.age i.male, breslow vce(cluster id)
        matrix bh = e(b)
        matrix vh = e(V)
        local nref = _N
    }
    frame drop Fr
    if `m' == 1 {
        matrix bref_a = ba / `M'
        matrix bref_t = bt / `M'
        matrix bref_h = bh / `M'
        matrix vref_h = vh / `M'
    }
    else {
        matrix bref_a = bref_a + ba / `M'
        matrix bref_t = bref_t + bt / `M'
        matrix bref_h = bref_h + bh / `M'
        matrix vref_h = vref_h + vh / `M'
    }
}
chk "pieces in the result = reference" `Np' `nref' 0.5
cmpw res flong "flong" `M' `Np'

* ===== 4. the weighted baseline table =================================================================================================================
di as res _n "=== 4. mi estimate: mean [pw=_ipw_ptw] (flong) ==="
frame res {
    capture noisily mi estimate, post: mean age biomarker [pw=_ipw_ptw] if adh == 1, over(arm)
    local rc_a = _rc
    matrix bmi_a = e(b)
    capture noisily mi estimate, post: mean age biomarker [pw=_ipw_ptwt] if adh == 1, over(arm)
    local rc_t = _rc
    matrix bmi_t = e(b)
}
chkrc "mi estimate: mean [pw=_ipw_ptw]" `rc_a' 0
chkrc "mi estimate: mean [pw=_ipw_ptwt]" `rc_t' 0
mata: st_numscalar("dtab_a", max(abs(st_matrix("bmi_a") - st_matrix("bref_a"))))
mata: st_numscalar("dtab_t", max(abs(st_matrix("bmi_t") - st_matrix("bref_t"))))
chk "mean table, ptw: max diff to the mean over m of per-m means" scalar(dtab_a) 0 1e-9
chk "mean table, ptwt: max diff" scalar(dtab_t) 0 1e-9
matrix list bmi_a, format(%9.4f)

* ===== 5. mi stset + mi estimate: stcox ==================================================================================================================
di as res _n "=== 5. mi stset [pw=_ipw_psw] + mi estimate: stcox (flong) ==="
frame res {
    capture noisily mi stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
    local rc_s = _rc
    capture noisily mi estimate, post: stcox i.arm c.age i.male, breslow vce(cluster id)
    local rc_h = _rc
    matrix bhr = e(b)
    matrix vhr = e(V)
    capture matrix vmi = e(W_mi)
    local rc_vw = _rc
}
chkrc "mi stset with the unregistered weight" `rc_s' 0
chkrc "mi estimate: stcox" `rc_h' 0
mata: st_numscalar("dhr", max(abs(st_matrix("bhr") - st_matrix("bref_h"))))
chk "pooled coefficients = mean over m of per-m stset/stcox" scalar(dhr) 0 1e-8
if `rc_vw' == 0 {
    mata: st_numscalar("dvw", max(abs(st_matrix("vmi") - st_matrix("vref_h"))))
    chk "within-imputation variance = mean of per-m V" scalar(dvw) 0 1e-9
}
else di as txt "  (e(W_mi) not available: within-imputation variance not compared)"
matrix list bhr, format(%9.5f)

* ===== 6. options: replace / change / notruncate / mistset / keep ===============================================================================
di as res _n "=== 6. options ==="
capture noisily ipwadh_mi, `opts' frame(res)
chkrc "second call without replace" `=_rc' 110
ipwadh_mi, `opts' frame(res, replace) mistset(psw) change
chktxt "current frame after change" "`c(frame)'" "res"
capture noisily mi estimate, post: stcox i.arm c.age i.male, breslow vce(cluster id)
chkrc "mi estimate: stcox after mistset(psw)" `=_rc' 0
matrix bhr2 = e(b)
mata: st_numscalar("dhr2", max(abs(st_matrix("bhr2") - st_matrix("bhr"))))
chk "mistset(psw) gives the same pooled coefficients" scalar(dhr2) 0 1e-10
frame change default

ipwadh_mi, `opts' frame(res, replace) notruncate
foreach v in _ipw_pswt _ipw_ptwt {
    capture frame res: confirm variable `v'
    chkrc "notruncate: no `v'" `=_rc' 111
}
frame res: quietly count if _mi_m == 1 & !missing(_ipw_ptw)
local nn1 = r(N)
mata: st_numscalar("nref1", sum(REFW1[., 3] :< .))
chk "notruncate: non-missing weights in m=1 = reference (the 16 zero-length pieces have no weight)" `nn1' scalar(nref1) 0.5

ipwadh_mi, `opts' frame(res, replace) keep(age male)
frame res {
    quietly ds
    local nv : word count `r(varlist)'
    di as txt "  variables kept with keep(age male): `r(varlist)'"
    capture confirm variable endreason
    chkrc "keep(): endreason dropped" `=_rc' 111
    capture confirm variable age
    chkrc "keep(): age kept" `=_rc' 0
}
cmpw res flong "keep()" `M' `Np'
capture frame drop res

* ===== 7. mistyle(wide) and flong / mlong SOURCES ========================================================================================================
di as res _n "=== 7a. mistyle(wide) (wide source) ==="
ipwadh_mi, `opts' frame(resw) mistyle(wide)
cmpw resw wide "wide result" `M' `Np'
frame resw {
    capture noisily mi estimate, post: mean age biomarker [pw=_ipw_ptw] if adh == 1, over(arm)
    local rc_w = _rc
    matrix bmi_w = e(b)
}
chkrc "wide result: mi estimate: mean" `rc_w' 0
mata: st_numscalar("dtab_w", max(abs(st_matrix("bmi_w") - st_matrix("bref_a"))))
chk "wide result: mean table = reference" scalar(dtab_w) 0 1e-9
capture frame drop resw

foreach st in flong mlong {
    di as res _n "=== 7b. source style `st' ==="
    mkmi `st' `M' `seed'
    capture noisily ipwadh_mi, `opts' frame(res_`st')
    chkrc "`st' source: ipwadh_mi" `=_rc' 0
    cmpw res_`st' flong "`st' source" `M' `Np'
    capture frame drop res_`st'
}
* NOTE: the reference weights come from the WIDE-source imputation (same seed; mi convert keeps the imputed values), so equal weights are expected.

di as res _n "test_ipwadh_mi_frame.do: $nchk checks, $nfail failed"
di as txt "Send me the whole log."
