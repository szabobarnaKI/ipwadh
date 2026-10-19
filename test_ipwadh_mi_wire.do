* test_ipwadh_mi_wire.do -- certification of the frame() branches of ipwadh.ado v0.2.1: mi data (hand-over to ipwadh_mi.ado v0.2.1) and ORDINARY data (copy + weightsonly in a new frame).  NOT YET RUN IN STATA.  Stata 16+ (run in Stata 19).
* Run from the folder containing ipwadh.ado (0.2.0), ipwadh_mi.ado (0.2.1), ipwadhc.ado, ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do, ipwadh_eval.ado and data/S8_N200/episodes.dta.
* Data: the first 100 persons of S8_N200, M = 2 imputations (small and quick: part of the quick set of cscript_ipwadh.do).
* ASSUMPTION (not yet checked): data/S8_N200 has the same variables as data/S8_N2000 (id t0 t1 adh arm inf dead ltfu age male tsi abx biomarker ninf), as used by test_ipwadh_mi_frame.do.
* What it certifies
*   1  error paths of the wiring: mi data without frame() (198), adjustse with mi data (198), mi-only options (mistyle, mistset, mikeep) and change without frame() on ordinary data (198 each); the user's mi data untouched
*   2  ipwadh, frame() == ipwadh_mi, frame() directly: row layout, weights, imputed covariates identical (default flong), for several option sets (notruncate, strata, mistyle(wide), MIKeep -> keep)
*   3  change switches to the result frame; MIStset(psw) gives the same pooled mi estimate: stcox as a manual mi stset + mi estimate on the result of ipwadh_mi
*   4  frame(, replace) and the refusal to overwrite (110); keepdata / weightsonly are accepted on mi data
*   5  frame() on ORDINARY data: the user's data untouched, result frame == ipwadh, weightsonly on a fresh copy (weights, stcox), refusals (110, 198), no frame left behind after a failed call, change
*   (the ordinary path WITHOUT frame() is certified by test_ipwadh_wrapper.do, which must still pass unchanged)
clear all
set more off
global nfail = 0
global nchk = 0
local M 2
local seed 20261011

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

* compare the listed variables of two result frames, row by row (same rows, same missing pattern, max abs diff)
capture program drop cmpfr
program define cmpfr
    args fa fb label vlist
    frame `fa': quietly count
    local na = r(N)
    frame `fb': quietly count
    local nb = r(N)
    chk "`label': rows of the two result frames" `na' `nb' 0.5
    if `na' == `nb' {
        frame `fa': mata: CWA = st_data(., "`vlist'")
        frame `fb': mata: CWB = st_data(., "`vlist'")
        mata: st_numscalar("dmiss", sum((CWA :>= .) :!= (CWB :>= .)))
        mata: st_numscalar("dval", max(abs(editmissing(CWA, 0) - editmissing(CWB, 0))))
        chk "`label': cells whose missing status differs" scalar(dmiss) 0 0.5
        chk "`label': max abs diff of `vlist'" scalar(dval) 0 1e-12
        capture mata: mata drop CWA CWB
    }
end

capture program drop mkmi
program define mkmi
    args style M seed
    use "data/S8_N200/episodes.dta", clear
    keep if id <= 100
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

* ordinary (non-mi) data for section 5: the raw S8_N200 has missing biomarker / abx (ipwadh refuses such rows, r(416)), so fill them with 0 (only the plumbing is tested here)
capture program drop mkord
program define mkord
    use "data/S8_N200/episodes.dta", clear
    quietly replace biomarker = 0 if missing(biomarker)
    quietly replace abx = 0 if missing(abx)
end

local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local opts "`base' failure(inf) `lists'"
local wv4 "_ipw_psw _ipw_pswt _ipw_ptw _ipw_ptwt"

mkmi wide `M' `seed'
quietly datasignature
local sig0 "`r(datasignature)'"
local N0 = _N

* ===== 1. error paths ==============================================================================================================================
di as res _n "=== 1. error paths ==="
capture noisily ipwadh, `opts'
chkrc "mi data without frame()" `=_rc' 198
capture noisily ipwadh, `opts' frame(w_err) adjustse
chkrc "adjustse with mi data" `=_rc' 198
capture noisily ipwadh, `opts' frame(w_err) mistyle(long)
chkrc "bad mistyle() reaches ipwadh_mi" `=_rc' 198
capture noisily ipwadh, `opts' frame(w_err) mistyle(wide) mistset(psw)
chkrc "mistset() with mistyle(wide)" `=_rc' 198
capture noisily ipwadh if biomarker > 0, `opts' frame(w_err)
chkrc "if on an imputed variable" `=_rc' 198
capture frame w_err: describe
global nchk = $nchk + 1
if _rc == 0 {
    global nfail = $nfail + 1
    di as err "  FAIL  a result frame w_err was left behind by a failed call"
}
else di as txt "  PASS  no frame left behind by failed calls"
quietly datasignature
chktxt "datasignature of the user's mi data unchanged" "`r(datasignature)'" "`sig0'"
chk "N of the user's data unchanged" _N `N0' 0.5

capture frame drop plain
frame create plain
frame plain: quietly use "data/S8_N200/episodes.dta", clear
capture noisily frame plain: ipwadh, `opts' mistyle(flong)
chkrc "mistyle() on ordinary data" `=_rc' 198
capture noisily frame plain: ipwadh, `opts' mistset(psw)
chkrc "mistset() on ordinary data" `=_rc' 198
capture noisily frame plain: ipwadh, `opts' mikeep(age)
chkrc "mikeep() on ordinary data" `=_rc' 198
capture noisily frame plain: ipwadh, `opts' change
chkrc "change without frame() on ordinary data" `=_rc' 198
capture frame drop plain

* ===== 2. ipwadh, frame() == ipwadh_mi, frame() ====================================================================================================
di as res _n "=== 2a. default (flong, truncate(1 99)) ==="
ipwadh, `opts' frame(wA)
chktxt "current frame after the call (no change option)" "`c(frame)'" "default"
ipwadh_mi, `opts' frame(wB)
frame wA {
    chktxt "wrapper result: mi style" "`_dta[_mi_style]'" "flong"
    chktxt "wrapper result: characteristic ipwadhmi_mistyle" "`_dta[ipwadhmi_mistyle]'" "flong"
}
cmpfr wA wB "default" "`wv4' biomarker abx id t0 t1"

di as res _n "=== 2b. notruncate + strata(male) ==="
ipwadh, `opts' strata(male) notruncate frame(wA, replace)
ipwadh_mi, `opts' strata(male) notruncate frame(wB, replace)
cmpfr wA wB "notruncate strata" "_ipw_psw _ipw_ptw biomarker abx id t0 t1"
capture frame wA: confirm variable _ipw_pswt
chkrc "notruncate: no _ipw_pswt in the wrapper result" `=_rc' 111

di as res _n "=== 2c. mistyle(wide) ==="
ipwadh, `opts' mistyle(wide) frame(wA, replace)
ipwadh_mi, `opts' mistyle(wide) frame(wB, replace)
frame wA {
    chktxt "wrapper result: mi style" "`_dta[_mi_style]'" "wide"
}
cmpfr wA wB "wide" "_1__ipw_psw _2__ipw_psw _1__ipw_ptwt _2__ipw_ptwt _1_biomarker _2_abx id t0 t1"

di as res _n "=== 2d. MIKeep(age male) is keep() of ipwadh_mi; verbose; keepdata / weightsonly accepted ==="
ipwadh, `opts' mikeep(age male) frame(wA, replace)
ipwadh_mi, `opts' keep(age male) frame(wB, replace)
cmpfr wA wB "mikeep" "`wv4' biomarker abx age male id t0 t1"
capture frame wA: confirm variable endreason
chkrc "mikeep(): endreason dropped (as keep() in ipwadh_mi)" `=_rc' 111
ipwadh, `opts' frame(wA, replace) keepdata
cmpfr wA wB "keepdata accepted" "`wv4' id t0 t1"
ipwadh, `opts' frame(wA, replace) weightsonly verbose
cmpfr wA wB "weightsonly verbose accepted" "`wv4' id t0 t1"

di as res _n "=== 2e. frame() exists: refused without replace (110), second call with replace works ==="
capture noisily ipwadh, `opts' frame(wA)
chkrc "existing frame without replace" `=_rc' 110
ipwadh, `opts' frame(wA, replace)
cmpfr wA wB "replace" "`wv4' id t0 t1"

* ===== 3. change and MIStset ========================================================================================================================
di as res _n "=== 3. change, MIStset(psw) and mi estimate: stcox ==="
ipwadh_mi, `opts' frame(wB, replace)
frame wB {
    quietly mi stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
    quietly mi estimate, post: stcox i.arm c.age i.male, breslow vce(cluster id)
    matrix bman = e(b)
}
ipwadh, `opts' frame(wA, replace) mistset(psw) change
chktxt "current frame after change" "`c(frame)'" "wA"
capture noisily mi estimate, post: stcox i.arm c.age i.male, breslow vce(cluster id)
chkrc "mi estimate: stcox after MIStset(psw) through ipwadh" `=_rc' 0
matrix bwr = e(b)
mata: st_numscalar("dhr", max(abs(st_matrix("bwr") - st_matrix("bman"))))
chk "pooled coefficients: wrapper+MIStset(psw) = manual mi stset on ipwadh_mi" scalar(dhr) 0 1e-10
matrix list bwr, format(%9.5f)
frame change default
quietly datasignature
chktxt "datasignature of the user's mi data still unchanged" "`r(datasignature)'" "`sig0'"

capture frame drop wA
capture frame drop wB

* ===== 4. flong source through the wrapper ============================================================================================================
di as res _n "=== 4. flong source ==="
mkmi flong `M' `seed'
capture noisily ipwadh, `opts' frame(wF)
chkrc "flong source through ipwadh" `=_rc' 0
mkmi wide `M' `seed'
ipwadh_mi, `opts' frame(wW)
cmpfr wF wW "flong source = wide source" "`wv4' biomarker abx id t0 t1"
capture frame drop wF
capture frame drop wW

* ===== 5. frame() on ORDINARY data ======================================================================================================================
di as res _n "=== 5. frame() on ordinary data ==="
mkord
quietly datasignature
local sigo "`r(datasignature)'"
local No = _N
ipwadh, `opts' frame(o1)
chktxt "current frame after the call (no change)" "`c(frame)'" "default"
quietly datasignature
chktxt "datasignature of the ordinary data unchanged" "`r(datasignature)'" "`sigo'"
chk "N of the ordinary data unchanged" _N `No' 0.5
capture confirm variable _ipw_psw
chkrc "no weights in the user's frame" `=_rc' 111
frame o1 {
    quietly count
    local Nsplit = r(N)
    di as txt "  pieces in the result frame: `Nsplit'"
}
global nchk = $nchk + 1
if `Nsplit' > `No' di as txt "  PASS  the result frame holds the split data (`Nsplit' > `No' rows)"
else {
    global nfail = $nfail + 1
    di as err "  FAIL  the result frame has `Nsplit' rows, not more than the unsplit data (`No')"
}

* the same pipeline directly (weightsonly on a fresh copy of the data)
capture frame drop o2
frame create o2
frame o2 {
    quietly mkord
    ipwadh, `opts' weightsonly
}
cmpfr o1 o2 "ordinary frame() == weightsonly" "_ipw_psw _ipw_pswt id t0 t1"
frame o1 {
    quietly stcox i.arm c.age i.male, breslow vce(cluster id)
    local bo1 = _b[1.arm]
}
frame o2 {
    quietly stcox i.arm c.age i.male, breslow vce(cluster id)
    local bo2 = _b[1.arm]
}
chk "stcox in the result frame (stset is there) = stcox after weightsonly" `bo1' `bo2' 1e-12
capture frame drop o2

* refusals and cleanup
capture noisily ipwadh, `opts' frame(o1)
chkrc "existing frame without replace" `=_rc' 110
capture noisily ipwadh, `opts' frame(default, replace)
chkrc "result frame = current frame" `=_rc' 198
capture noisily ipwadh, `opts' frame(o3) adjustse
chkrc "adjustse with frame()" `=_rc' 198
capture noisily ipwadh, `opts' frame(o3) mistyle(wide)
chkrc "mistyle() with ordinary data" `=_rc' 198
capture noisily ipwadh, `opts' frame(o3) truncate(1 99) notruncate
chkrc "truncate + notruncate through frame()" `=_rc' 198
capture noisily frame o1: ipwadh, `opts' frame(o3)
chkrc "frame() from already split data" `=_rc' 110
capture frame o3: describe
global nchk = $nchk + 1
if _rc == 0 {
    global nfail = $nfail + 1
    di as err "  FAIL  a result frame o3 was left behind by a failed call"
}
else di as txt "  PASS  no frame left behind by failed calls"

* replace, change, notruncate
ipwadh, `opts' frame(o1, replace) notruncate change
chktxt "current frame after change" "`c(frame)'" "o1"
capture confirm variable _ipw_pswt
chkrc "notruncate: no _ipw_pswt in the result frame" `=_rc' 111
capture confirm variable _ipw_psw
chkrc "_ipw_psw in the result frame" `=_rc' 0
frame change default
capture frame drop o1
quietly datasignature
chktxt "datasignature of the ordinary data still unchanged" "`r(datasignature)'" "`sigo'"

di as res _n "test_ipwadh_mi_wire.do: $nchk checks, $nfail failed"
di as txt "Send me the whole log."
