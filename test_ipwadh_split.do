* test_ipwadh_split.do -- run in Stata 19 from the folder containing ipwadh_split.ado, ipwadh_mata.do and data/
* NOT YET RUN. Expected numbers come from the Python reference (ipwadh_simdata.py), dataset S3_N200 (seed 1):
*   exact default 77,718 rows; allrows 125,641; |C| = 3,425 (max bandpoints 1,712); bandpoints(10) 42,062; bandpoints(100) 10,918;
*   cutall 430,066; strata(male) 42,001; strata(arm male) 24,165; strata(arm male)+bandpoints(10) 21,264;
*   auto: maxrows(100000) -> m=1, maxrows(60000) -> m=5 (57,592 rows), maxrows(40000) -> m=11 (39,396), maxrows(20000) -> m=32 (19,783),
*   maxrows(8000) -> m=328 (7,987); maxrows(5000) infeasible (rc 950).
*   force: overrides maxgb()/maxrows()/set max_memory, never the observation limit (tested with the testing aid maxobs(50000)).
clear all
set more off
local d "data/S3_N200"
local base "id(id) time0(t0) time1(t1) adherence(adh) failure(inf)"

* ---- 1. pre-check only: row counts ------------------------------------------------------------
use "`d'/episodes.dta", clear
local N0 = _N
ipwadh_split, `base' check
assert r(N_after) == 77718 & r(exact) == 1 & r(m) == 1
assert r(n_cutpoints) == 3425 & r(max_bandpoints) == 1712 & r(N_exact) == 77718
ipwadh_split, `base' allrows check
assert r(N_after) == 125641
ipwadh_split, `base' bandpoints(10) check
assert r(N_after) == 42062 & r(exact) == 0 & r(N_exact) == 77718
ipwadh_split, `base' bandpoints(100) check
assert r(N_after) == 10918
ipwadh_split, `base' cutall check
assert r(N_after) == 430066
ipwadh_split, `base' strata(male) check
assert r(N_after) == 42001
ipwadh_split, `base' strata(arm male) check
assert r(N_after) == 24165
ipwadh_split, `base' strata(arm male) bandpoints(10) check
assert r(N_after) == 21264
assert _N == `N0'                       // check must not change the data
di as res "check: row counts match Python reference"

* ---- 2. bandpoints(auto) -----------------------------------------------------------------------------
foreach c in "100000 1 77718" "60000 5 57592" "40000 11 39396" "20000 32 19783" "8000 328 7987" {
    tokenize `c'
    ipwadh_split, `base' bandpoints(auto) maxrows(`1') check
    assert r(m) == `2' & r(N_after) == `3'
    assert r(N_after) <= `1'
}
capture ipwadh_split, `base' bandpoints(auto) maxrows(5000) check
assert _rc == 950
assert _N == `N0'
di as res "auto: m and row counts match Python reference; infeasible budget refused"

* ---- 3. option validation ------------------------------------------------------------------------------
capture ipwadh_split, `base' bandpoints(1713) check
assert _rc == 198
capture ipwadh_split, `base' bandpoints(0) check
assert _rc == 198
capture ipwadh_split, `base' bandpoints(abc) check
assert _rc == 198
ipwadh_split, `base' bandpoints(1712) check                  // allowed: exactly 50%
assert _N == `N0'
di as res "bandpoints() validation ok"

* ---- 4. guard: exact split too big for maxgb(): refuse, data untouched (fixed bandpoints never bands silently) ---------
capture noisily ipwadh_split, `base' maxgb(0.0001)
assert _rc == 950
assert _N == `N0'
di as res "guard: refused without force"

* ---- 4b. force overrides every MEMORY limit but never the observation limit ---------------------------------------------
ipwadh_split, `base' maxgb(0.0001) force check
assert r(ok) == 1 & r(forced) == 1 & r(N_after) == 77718 & r(m) == 1
ipwadh_split, `base' maxrows(1000) force check                   // maxrows() is also a memory-type limit
assert r(ok) == 1 & r(forced) == 1
ipwadh_split, `base' bandpoints(auto) maxgb(0.0001) force check  // force removes the memory budget: auto then gives m=1
assert r(m) == 1 & r(N_after) == 77718
ipwadh_split, `base' maxobs(50000) force check                   // observation limit: not overridable
assert r(ok) == 0
capture ipwadh_split, `base' maxobs(50000) force
assert _rc == 900
assert _N == `N0'
ipwadh_split, `base' maxobs(50000) bandpoints(auto) check        // auto respects the observation limit by banding
assert r(ok) == 1 & r(m) == 8 & r(N_after) == 47852      // Python reference: m=8, 47,852 rows
* forced split really runs and gives the exact result
use "`d'/episodes.dta", clear
ipwadh_split, `base' maxgb(0.0001) force
assert _N == 77718
di as res "force: memory limits overridden, observation limit respected"

* ---- 5. split and compare piece by piece with the Python reference ------------------------------------------------------
local cases `""split_adherent|" "split_adherent_m10|bandpoints(10)" "split_adherent_strata_male|strata(male)" "split_adherent_strata_armmale_m10|strata(arm male) bandpoints(10)""'
foreach c of local cases {
    gettoken ref opts : c, parse("|")
    local opts = subinstr("`opts'", "|", "", 1)
    use "`d'/episodes.dta", clear
    ipwadh_split, `base' `opts'
    gen double dur = t1 - t0
    quietly summarize dur
    local pt = r(sum)
    quietly summarize inf
    local ninf = r(sum)
    keep id _ipw_orig _ipw_piece t0 t1 inf
    rename (t0 t1 inf) (s_t0 s_t1 s_inf)
    tempfile mine
    save `mine', replace
    use "`d'/`ref'.dta", clear
    bysort orig_row: gen long _ipw_piece = _n
    rename orig_row _ipw_orig
    merge 1:1 _ipw_orig _ipw_piece using `mine', assert(match) nogenerate
    assert reldif(t0, s_t0) < 1e-12 & reldif(t1, s_t1) < 1e-12
    assert inf == s_inf
    use "`d'/episodes.dta", clear
    gen double dur = t1 - t0
    quietly summarize dur
    assert reldif(r(sum), `pt') < 1e-12           // person-time conserved
    quietly summarize inf
    assert r(sum) == `ninf'                       // infections conserved
    di as res "`ref': pieces identical to Python reference; person-time and infections conserved"
}
di as res "test_ipwadh_split.do: all assertions passed"
