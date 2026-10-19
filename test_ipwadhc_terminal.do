* test_ipwadhc_terminal.do -- run in Stata 19 from the folder containing ipwadhc.ado (v0.3.0), ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do,
* ipwadh_eval.ado (v0.2.0), expected_ipwadhc_terminal_S3_N200.csv and data/S3_N200/episodes.dta, data/S3_N2000/episodes.dta.   NOT YET RUN.
* Tests ipwadhc 0.3.0: option ltfu() (second censoring cause, IPCW) and terminal(censor|none).  Expected values: independent Python reference
* (ref_core.py / ref_terminal_expected.py / ref_terminal.py, Session 5a; Poisson fits converged to 1e-12, same models as the Stata call).
*   Sections 1-4  N=200 row by row (6,754 usable rows): death+ltfu, ltfu only, death only, terminal(none), no censoring; parts identity; returns
*   Section 5     N=2000 summaries of the row-end weights
*   Section 6     error handling (198 / 459)
*   Section 7     N=2000 full pipeline: death+ltfu IPCW (exact split, evalat end, truncated) and the composite analysis (failure = infection or death)
* Checks do NOT stop at the first difference: PASS/FAIL per check, totals at the end. Send me the whole log.  Section 7 needs about 2 GB and 5-8 minutes.
clear all
set more off
global nfail = 0
global nchk = 0

capture program drop chk
program define chk
    args label value expected tol mode
    if "`mode'" == "" local mode "abs"
    if "`mode'" == "abs"  local dif = abs(`value' - `expected')
    if "`mode'" == "rel"  local dif = reldif(`value', `expected')
    if "`mode'" == "rel2" local dif = abs(`value' - `expected') / abs(`expected')
    global nchk = $nchk + 1
    if `dif' < `tol' {
        di as txt "  PASS  `label': " %16.10f `value' "  (python " %16.10f `expected' ")"
    }
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label': " %16.10f `value' "  (python " %16.10f `expected' ", `mode' diff " %10.3e `dif' ")"
    }
end

* row-level comparison of the weights in memory with the Python columns merged from the expected file
capture program drop cmp_rows
program define cmp_rows
    args label swc_py sw_py
    preserve
    merge 1:1 id rowno using "expected_terminal_tmp.dta", keep(master match) generate(mg)
    quietly count if mg == 3
    chk "`label': rows matched with the Python reference" r(N) 6754 0.5
    quietly count if (mg == 1) != (t1 == t0)
    chk "`label': unmatched rows are exactly the zero-length rows" r(N) 0 0.5
    quietly gen double rd_a = abs(_ipw_swa - x_swa) / abs(x_swa) if mg == 3
    quietly gen double rd_c = abs(_ipw_swc - `swc_py') / abs(`swc_py') if mg == 3
    quietly gen double rd_s = abs(_ipw_sw - `sw_py') / abs(`sw_py') if mg == 3
    foreach v in a c s {
        quietly summarize rd_`v'
        local mx = r(max)
        chk "`label': row-level max relative difference, `v'" `mx' 0 1e-5
    }
    restore
end

* expected row-level weights from the Python reference
import delimited using "expected_ipwadhc_terminal_S3_N200.csv", clear
foreach v in swa swc_d swc_l swc_dl sw_d sw_l sw_dl {
    rename `v' x_`v'
}
assert _N == 6763                      // the file has all rows; the 9 zero-length rows have missing weights
drop if missing(x_swa)
assert _N == 6754
save "expected_terminal_tmp.dta", replace

local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm)"
timer clear
timer on 1

* ===== 1. N=200: death + loss to follow-up (terminal(censor), the default) ===================================================================
di as res _n "=== 1. death + ltfu, default terminal(censor), parts ==="
use "data/S3_N200/episodes.dta", clear
quietly count if ltfu == 1
chk "N200 ltfu events in the data" r(N) 23 0.5
ipwadhc, `base' dead(dead) ltfu(ltfu) `lists' parts replace
assert "`r(terminal)'" == "censor"
assert "`r(censoring)'" == "death ltfu"
assert r(n_cens_models) == 2
chk "N200 usable rows" r(N) 6754 0.5
cmp_rows "N200 death+ltfu" x_swc_dl x_sw_dl
* parts identity: swc = exp(Ec0 - dC * row length), swa = exp(E0 - dA * len + J0)
quietly gen double chk_c = abs(_ipw_swc - exp(_ipw_Ec0 - _ipw_dC * (t1 - t0))) if t1 > t0
quietly summarize chk_c
chk "N200 parts identity swc = exp(Ec0 - dC len), max abs diff" r(max) 0 1e-12
quietly gen double chk_a = abs(_ipw_swa - exp(_ipw_E0 - _ipw_dA * (t1 - t0) + _ipw_J0)) if t1 > t0
quietly summarize chk_a
chk "N200 parts identity swa, max abs diff" r(max) 0 1e-12
quietly summarize _ipw_swc if t1 > t0
chk "N200 swc (death+ltfu) mean" r(mean) 1.00332979 1e-5 rel2
chk "N200 swc (death+ltfu) max" r(max) 1.16232050 1e-5 rel2
chk "N200 swc (death+ltfu) min" r(min) 0.93359444 1e-5 rel2

* ===== 2. N=200: loss to follow-up only (no dead()) ============================================================================================
di as res _n "=== 2. ltfu only ==="
use "data/S3_N200/episodes.dta", clear
ipwadhc, `base' ltfu(ltfu) `lists' replace
assert "`r(censoring)'" == "ltfu"
assert r(n_cens_models) == 1
cmp_rows "N200 ltfu only" x_swc_l x_sw_l

* ===== 3. N=200: death only (as in 0.2.3), and terminal(none) with dead() and ltfu() =========================================================
di as res _n "=== 3. death only; terminal(none) ==="
use "data/S3_N200/episodes.dta", clear
ipwadhc, `base' dead(dead) `lists' replace
assert "`r(censoring)'" == "death"
cmp_rows "N200 death only" x_swc_d x_sw_d
use "data/S3_N200/episodes.dta", clear
ipwadhc, `base' dead(dead) ltfu(ltfu) terminal(none) `lists' replace
assert "`r(terminal)'" == "none"
assert "`r(censoring)'" == "ltfu"
cmp_rows "N200 terminal(none) + ltfu (= ltfu only)" x_swc_l x_sw_l
use "data/S3_N200/episodes.dta", clear
ipwadhc, `base' dead(dead) terminal(none) `lists' replace
assert r(n_cens_models) == 0
quietly summarize _ipw_swc
chk "N200 terminal(none), no ltfu: swc min" r(min) 1 1e-12
chk "N200 terminal(none), no ltfu: swc max" r(max) 1 1e-12

* ===== 4. N=200: no censoring model at all ======================================================================================================
di as res _n "=== 4. neither dead() nor ltfu() ==="
use "data/S3_N200/episodes.dta", clear
ipwadhc, `base' `lists' parts replace
assert r(n_cens_models) == 0
quietly summarize _ipw_swc
chk "N200 no censoring models: swc min" r(min) 1 1e-12
chk "N200 no censoring models: swc max" r(max) 1 1e-12
quietly summarize _ipw_dC
local mxabs = max(abs(r(min)), abs(r(max)))
chk "N200 no censoring models: dC max abs" `mxabs' 0 1e-12

* ===== 5. N=2000 summaries of the row-end weights (Python expected) ============================================================================
di as res _n "=== 5. N=2000 summaries ==="
use "data/S3_N2000/episodes.dta", clear
ipwadhc, `base' dead(dead) ltfu(ltfu) `lists' replace
chk "N2000 usable rows" r(N) 69572 0.5
foreach spec in "swc 1.00227663 0.01984890 0.82894675 1.15911569" "sw 0.98691495 1.35597159 0.00005519 41.16740191" {
    local v : word 1 of `spec'
    local m : word 2 of `spec'
    local s : word 3 of `spec'
    local mn : word 4 of `spec'
    local mx : word 5 of `spec'
    quietly summarize _ipw_`v' if t1 > t0
    chk "N2000 death+ltfu `v' mean" r(mean) `m' 1e-4 rel2
    chk "N2000 death+ltfu `v' sd" r(sd) `s' 1e-4 rel2
    chk "N2000 death+ltfu `v' min" r(min) `mn' 1e-4 rel2
    chk "N2000 death+ltfu `v' max" r(max) `mx' 1e-4 rel2
}
use "data/S3_N2000/episodes.dta", clear
ipwadhc, `base' ltfu(ltfu) `lists' replace
quietly summarize _ipw_swc if t1 > t0
chk "N2000 ltfu only swc mean" r(mean) 1.00085339 1e-4 rel2
chk "N2000 ltfu only swc sd" r(sd) 0.00431187 1e-4 rel2
chk "N2000 ltfu only swc max" r(max) 1.04119018 1e-4 rel2

* ===== 6. errors ================================================================================================================================
di as res _n "=== 6. error handling ==="
use "data/S3_N200/episodes.dta", clear
capture noisily ipwadhc, `base' dead(dead) terminal(foo) `lists' replace
chk "terminal(foo) is an error, rc" _rc 198 0.5
capture noisily ipwadhc, `base' dead(dead) ltfu(dead) `lists' replace
chk "dead() and ltfu() the same variable, rc" _rc 198 0.5
* both 1 on one row
gen byte l2 = ltfu
quietly summarize id if dead == 1, meanonly
quietly replace l2 = 1 if dead == 1 & id == r(min)
capture noisily ipwadhc, `base' dead(dead) ltfu(l2) `lists' replace
chk "dead and ltfu both 1 on a row, rc" _rc 459 0.5
drop l2
* ltfu = 1 on a row that is not the last usable row of the person
bysort id (t0 t1): gen byte cc = (_n == 2 & _N > 3 & t1 > t0)
quietly summarize id if cc == 1, meanonly
gen byte l3 = (cc == 1 & id == r(min))
capture noisily ipwadhc, `base' ltfu(l3) `lists' replace
chk "ltfu = 1 before the last row, rc" _rc 459 0.5
drop cc l3
gen byte l4 = 2 in 1
capture noisily ipwadhc, `base' ltfu(l4) `lists' replace
chk "ltfu not coded 0/1, rc" _rc 198 0.5
drop l4
* no ltfu events at all: the model cannot be fitted
gen byte l5 = 0
capture noisily ipwadhc, `base' ltfu(l5) `lists' replace
chk "ltfu() without any event, rc" _rc 2000 0.5

* ===== 7. N=2000 pipeline ========================================================================================================================
* Python (ref_terminal.py, exact adherent-only split, evalat end, Stata cluster SE):
*   death + ltfu IPCW : coef -0.391823  SE 0.053725 ; truncated p1/p99 (person-time) coef -0.390331  SE 0.050138
*   composite (failure = infection or death), adherence weights only: 7,794,187 rows, coef -0.372770  SE 0.052887 ; truncated -0.376275 / 0.049023
*   composite, unweighted: -0.370471 / 0.040822 ; infection only, unweighted -0.379687 / 0.041766
di as res _n "=== 7a. N2000 pipeline, death + ltfu IPCW ==="
use "data/S3_N2000/episodes.dta", clear
ipwadhc, `base' dead(dead) ltfu(ltfu) `lists' parts replace
ipwadh_split, id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu)
chk "exact split rows" _N 7524605 0.5
ipwadh_eval, time0(t0) time1(t1) adherence(adh) truncate(1 99) replace
stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
stcox i.arm, breslow vce(cluster id)
chk "death+ltfu IPCW coef" _b[1.arm] -0.391823 2e-4
chk "death+ltfu IPCW SE" _se[1.arm] 0.053725 2e-4
stset t1 [pw=_ipw_pswt] if adh == 1, time0(t0) failure(inf)
stcox i.arm, breslow vce(cluster id)
chk "death+ltfu IPCW truncated coef" _b[1.arm] -0.390331 2e-4
chk "death+ltfu IPCW truncated SE" _se[1.arm] 0.050138 2e-4

di as res _n "=== 7b. N2000 composite analysis (infection or death), terminal(none) ==="
use "data/S3_N2000/episodes.dta", clear
gen byte evt = (inf == 1 | dead == 1)
ipwadhc, `base' dead(dead) terminal(none) `lists' parts replace
ipwadh_split, id(id) time0(t0) time1(t1) adherence(adh) failure(evt) lastonly(inf dead ltfu)
chk "composite exact split rows" _N 7794187 0.5
ipwadh_eval, time0(t0) time1(t1) adherence(adh) truncate(1 99) replace
stset t1 if adh == 1, time0(t0) failure(evt)
stcox i.arm, breslow vce(cluster id)
chk "composite unweighted coef" _b[1.arm] -0.370471 2e-4
chk "composite unweighted SE" _se[1.arm] 0.040822 2e-4
stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(evt)
stcox i.arm, breslow vce(cluster id)
chk "composite adherence-weighted coef" _b[1.arm] -0.372770 2e-4
chk "composite adherence-weighted SE" _se[1.arm] 0.052887 2e-4
stset t1 [pw=_ipw_pswt] if adh == 1, time0(t0) failure(evt)
stcox i.arm, breslow vce(cluster id)
chk "composite adherence-weighted truncated coef" _b[1.arm] -0.376275 2e-4
chk "composite adherence-weighted truncated SE" _se[1.arm] 0.049023 2e-4

erase "expected_terminal_tmp.dta"
timer off 1
timer list 1
di as res _n "test_ipwadhc_terminal.do: $nchk checks, $nfail failed"
