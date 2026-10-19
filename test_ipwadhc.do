* test_ipwadhc.do -- run in Stata 19 from the folder containing ipwadhc.ado, ipwadh_covlists.ado and data/.   NOT YET RUN.
* Checks that need no external reference: zero-length rows, invariance (dropping ignored rows, row order, time unit), consistency of the
* parts() pieces with the row-end weight, input errors. Data: S3_N200 (6,763 rows; has zero-length rows and a rowno key).
clear all
set more off
local d "data/S3_N200"
local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"

use "`d'/episodes.dta", clear
quietly count if t1 == t0
local nz = r(N)
di as txt "zero-length rows in the data: `nz'"
assert `nz' > 0

* ---- 1. basic run: zero-length rows have missing weights, all others valid and positive -----------------------------------
ipwadhc, `base' `lists' parts replace
local rn = r(N)
local rz = r(n_zero)
local wv "`r(weightvar)'"
assert `rz' == `nz'
assert `rn' == _N - `nz'
assert missing(_ipw_swa) & missing(_ipw_sw) if t1 == t0
assert !missing(_ipw_swa) & !missing(_ipw_swc) & !missing(_ipw_sw) if t1 > t0
assert _ipw_swa > 0 & _ipw_swc > 0 & _ipw_sw > 0 if t1 > t0
assert "`wv'" == "_ipw_sw"
di as res "1. zero-length rows ignored, other weights valid: ok"

* ---- 2. parts reproduce the row-end weight:  swa = exp(E0 - dA*len + J0),  swc = exp(Ec0 - dC*len) ---------------------------
gen double len_ = t1 - t0
assert reldif(_ipw_swa, exp(_ipw_E0 - _ipw_dA * len_ + _ipw_J0)) < 1e-10 if t1 > t0
assert reldif(_ipw_swc, exp(_ipw_Ec0 - _ipw_dC * len_)) < 1e-10 if t1 > t0
* first usable row of each person: no exposure before it, no jump before it
bysort id (t0 t1 rowno): gen long k_ = sum(t1 > t0) if t1 > t0
assert _ipw_E0 == 0 & _ipw_J0 == 0 & _ipw_Ec0 == 0 if k_ == 1
drop len_ k_
di as res "2. parts consistent with the row-end weights: ok"

preserve
keep id rowno _ipw_swa _ipw_swc _ipw_sw
rename (_ipw_swa _ipw_swc _ipw_sw) (a1 c1 s1)
tempfile r1
save `r1'
restore

* ---- 3. dropping the ignored (zero-length) rows changes nothing -------------------------------------------------------------
preserve
drop if t1 == t0
ipwadhc, `base' `lists' replace
merge 1:1 id rowno using `r1', assert(match using) keep(match) nogenerate
assert reldif(_ipw_swa, a1) < 1e-10 & reldif(_ipw_swc, c1) < 1e-10 & reldif(_ipw_sw, s1) < 1e-10
restore
di as res "3. result identical after dropping zero-length rows: ok"

* ---- 4. row order does not matter ---------------------------------------------------------------------------------------------
preserve
set seed 20261009
gen double u_ = runiform()
sort u_
drop u_
ipwadhc, `base' `lists' replace
merge 1:1 id rowno using `r1', assert(match) nogenerate
assert reldif(_ipw_swa, a1) < 1e-10 if t1 > t0
assert reldif(_ipw_sw, s1) < 1e-10 if t1 > t0
assert missing(_ipw_sw) if t1 == t0
restore
di as res "4. row order invariance: ok"

* ---- 5. time unit does not matter (days -> years; rates scale, exposure scales back) -----------------------------------------
preserve
replace t0 = t0 / 365.25
replace t1 = t1 / 365.25
ipwadhc, `base' `lists' replace
merge 1:1 id rowno using `r1', assert(match) nogenerate
assert reldif(_ipw_swa, a1) < 1e-4 if t1 > t0
assert reldif(_ipw_sw, s1) < 1e-4 if t1 > t0
restore
di as res "5. time-unit invariance: ok"

* ---- 6. truncation and byarm run; truncated weights stay within the percentiles ---------------------------------------------
preserve
ipwadhc, `base' `lists' truncate(1 99) replace
local wv "`r(weightvar)'"
assert "`wv'" == "_ipw_swt"
quietly _pctile _ipw_sw if t1 > t0, p(1 99)
local lo = r(r1)
local hi = r(r2)
assert reldif(_ipw_swt, min(max(_ipw_sw, `lo'), `hi')) < 1e-12 if t1 > t0
assert missing(_ipw_swt) if t1 == t0
ipwadhc, `base' `lists' byarm replace
assert _ipw_sw > 0 if t1 > t0
restore
di as res "6. truncate and byarm run: ok"

* ---- 7. no death model: censoring weight is 1 -------------------------------------------------------------------------------
preserve
ipwadhc, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) `lists' parts replace
assert _ipw_swc == 1 & _ipw_Ec0 == 0 & _ipw_dC == 0 if t1 > t0
restore
di as res "7. without dead(): swc == 1: ok"

* ---- 8. input errors --------------------------------------------------------------------------------------------------------
preserve
replace adh = 2 in 1
capture ipwadhc, `base' `lists' replace
assert _rc == 198
restore

preserve
replace t1 = t0 - 1 in 1                           // negative length on row 1 (its t0 is not missing)
capture ipwadhc, `base' `lists' replace
assert _rc == 459
restore

capture ipwadhc, `base' stcoxvars(c.age) numvars(c.age) replace      // plain-language and advanced options mixed
assert _rc == 198
capture ipwadhc, `base' `lists'                                      // _ipw_* variables from section 1 exist, no -replace-
assert _rc == 110
di as res "8. input errors: ok"

* ---- 9. export the weights for the row-by-row comparison with the Python reference (ref_ipwadhc.py) ------------------------------
ipwadhc, `base' `lists' replace
export delimited id rowno _ipw_swa _ipw_swc _ipw_sw using "stata_weights_S3_N200.csv", replace
di as res "9. exported stata_weights_S3_N200.csv (send it back for: python3 ref_ipwadhc.py data/S3_N200/episodes.csv stata_weights_S3_N200.csv)"

di as res "test_ipwadhc.do: all assertions passed"
