* test_ipwadh_eval.do -- run in Stata 19 from the folder containing ipwadhc.ado, ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do,
* ipwadh_eval.ado and data/.   NOT YET RUN.   Data: S3_N200 (6,763 episode rows; 9 zero-length rows; 6,754 usable).
* Expected numbers come from ref_piece_weights.py (independent Python reference; same weights up to the Poisson optimizer tolerance ~6e-7):
*   exact split (m=1) 77,718 pieces, evalat(end):  coef -0.4711600819  SE(Stata vce(cluster)) 0.1821656452  (G = 192 clusters)
*   bandpoints(10) 42,062 pieces, evalat(mid):     coef -0.4722568068  SE 0.1820087318 ; evalat(end): -0.4713298139 / 0.1818826235
clear all
set more off
local d "data/S3_N200"
local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local sp "id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu)"

* =========================== A. exact split (bandpoints(1)) ===========================
use "`d'/episodes.dta", clear
quietly count if t1 == t0
local nz = r(N)
assert `nz' == 9
ipwadhc, `base' `lists' parts replace
ipwadh_split, `sp'
local m_ = r(m)
assert `m_' == 1
assert _N == 77718

* ---- 1. default evalat(end): zero-length rows have missing weights, all others valid and positive ---------------------------------
ipwadh_eval, time0(t0) time1(t1) adherence(adh) replace
local nmiss_ = r(n_missing)
local wv_ = "`r(weightvar)'"
local ev_ = "`r(evalat)'"
assert `nmiss_' == `nz'
assert "`wv_'" == "_ipw_psw"
assert "`ev_'" == "end"
assert missing(_ipw_psw) if t1 == t0
assert !missing(_ipw_pswa) & !missing(_ipw_pswc) & !missing(_ipw_psw) if t1 > t0
assert _ipw_pswa > 0 & _ipw_pswc > 0 & _ipw_psw > 0 if t1 > t0
assert reldif(_ipw_psw, _ipw_pswa * _ipw_pswc) < 1e-14 if t1 > t0
di as res "1. evalat(end): missing only for zero-length rows, others positive: ok"

* ---- 2. last piece of every row at evalat(end) reproduces the row-end weight of ipwadhc (_ipw_sw = swa*swc) ----------------------
bysort _ipw_orig (_ipw_piece): gen byte last_ = (_n == _N)
assert reldif(_ipw_psw, _ipw_sw) < 1e-10 if last_ == 1 & t1 > t0
assert reldif(_ipw_pswa, _ipw_swa) < 1e-10 & reldif(_ipw_pswc, _ipw_swc) < 1e-10 if last_ == 1 & t1 > t0
di as res "2. last piece, evalat(end) == row-end weight from ipwadhc: ok"

* ---- 3. log-linearity: start, mid, end ------------------------------------------------------------------------------------------
ipwadh_eval, time0(t0) time1(t1) evalat(start) generate(_st) replace
ipwadh_eval, time0(t0) time1(t1) evalat(mid)   generate(_md) replace
* start of the first piece = weight at the start of the row
assert reldif(_st_psw, exp(_ipw_E0 + _ipw_J0 + _ipw_Ec0)) < 1e-10 if _ipw_piece == 1 & t1 > t0
* continuity: start of piece j+1 = end of piece j
bysort _ipw_orig (_ipw_piece): assert reldif(_st_psw, _ipw_psw[_n-1]) < 1e-10 if _n > 1 & !missing(_st_psw)
* log W(mid) = (log W(start) + log W(end))/2
assert reldif(_md_psw, sqrt(_st_psw * _ipw_psw)) < 1e-10 if t1 > t0
di as res "3. start / mid / end: continuity and log-linearity: ok"

* ---- 4. export pieces for the Python comparison, then the on-treatment Cox fit --------------------------------------------------
preserve
keep id rowno _ipw_piece _ipw_pswa _ipw_pswc _ipw_psw
export delimited using "stata_pw_m1_end.csv", replace
restore

stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)     // NO id(): Stata requires weights constant within id()
stcox i.arm, breslow vce(cluster id)
di "m=1 end:  coef " %12.8f _b[1.arm] "  SE " %12.8f _se[1.arm]
assert abs(_b[1.arm]  - (-0.4711600819)) < 2e-5
assert abs(_se[1.arm] - 0.1821656452)   < 2e-5
di as res "4. exact split, evalat(end): Cox coef and SE match the Python reference: ok"

* ---- 5. row order does not matter -------------------------------------------------------------------------------------------------
preserve
keep _ipw_orig _ipw_piece _ipw_psw
rename _ipw_psw ref_
tempfile r1
save `r1'
restore
preserve
set seed 20261009
gen double u_ = runiform()
sort u_
drop u_
ipwadh_eval, time0(t0) time1(t1) replace
merge 1:1 _ipw_orig _ipw_piece using `r1', assert(match) nogenerate
assert reldif(_ipw_psw, ref_) < 1e-12 if t1 > t0
assert missing(_ipw_psw) if t1 == t0
restore
di as res "5. row order invariance: ok"

* ---- 6. input errors -------------------------------------------------------------------------------------------------------------
capture ipwadh_eval, time0(t0) time1(t1) evalat(middle) replace
assert _rc == 198
capture ipwadh_eval, time0(t0) time1(t1)                      // _ipw_psw etc. exist, no -replace-
assert _rc == 110
preserve
quietly levelsof _ipw_orig if _ipw_piece == 3, local(o3)
local o : word 1 of `o3'
drop if _ipw_orig == `o' & _ipw_piece == 2                    // damaged split
capture ipwadh_eval, time0(t0) time1(t1) replace
assert _rc == 459
restore
preserve
use "`d'/episodes.dta", clear                                 // unsplit data, no parts, no _ipw_orig
capture ipwadh_eval, time0(t0) time1(t1)
assert _rc == 111
restore
di as res "6. input errors (198, 110, 459, 111): ok"

* =========================== B. banded split, bandpoints(10) ===========================
use "`d'/episodes.dta", clear
ipwadhc, `base' `lists' parts replace
ipwadh_split, `sp' bandpoints(10)
local m_ = r(m)
assert `m_' == 10
assert _N == 42062
ipwadh_eval, time0(t0) time1(t1) evalat(mid) adherence(adh) replace
preserve
keep id rowno _ipw_piece _ipw_pswa _ipw_pswc _ipw_psw
export delimited using "stata_pw_m10_mid.csv", replace
restore
stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)     // NO id(): Stata requires weights constant within id()
stcox i.arm, breslow vce(cluster id)
di "m=10 mid: coef " %12.8f _b[1.arm] "  SE " %12.8f _se[1.arm]
assert abs(_b[1.arm]  - (-0.4722568068)) < 2e-5
assert abs(_se[1.arm] - 0.1820087318)   < 2e-5

ipwadh_eval, time0(t0) time1(t1) evalat(end) replace
preserve
keep id rowno _ipw_piece _ipw_pswa _ipw_pswc _ipw_psw
export delimited using "stata_pw_m10_end.csv", replace
restore
stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)     // NO id(): Stata requires weights constant within id()
stcox i.arm, breslow vce(cluster id)
di "m=10 end: coef " %12.8f _b[1.arm] "  SE " %12.8f _se[1.arm]
assert abs(_b[1.arm]  - (-0.4713298139)) < 2e-5
assert abs(_se[1.arm] - 0.1818826235)   < 2e-5
di as res "7. bandpoints(10): evalat(mid) and evalat(end) match the Python reference: ok"

di as res "test_ipwadh_eval.do: all assertions passed"
di as txt "Send back (or run yourself) for the piece-by-piece comparison:"
di as txt "  python3 ref_piece_weights.py data/S3_N200/episodes.csv stata_pw_m1_end.csv:1:end stata_pw_m10_mid.csv:10:mid stata_pw_m10_end.csv:10:end"
