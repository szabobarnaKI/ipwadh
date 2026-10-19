* cscript_ipwadh.do -- certification script for the ipwadh command suite
* ---------------------------------------------------------------------------------------------------------------------------------------
* WHAT IT DOES   Runs every test do-file of the suite in a clean Stata state, tallies pass / fail, writes the full log to cscript_ipwadh.log
*                and prints a one-page summary.  Stops with an error (r(9)) at the end if anything failed, so a batch run flags it.
* COMMANDS TESTED  ipwadhc (+ ipwadh_covlists), ipwadh_split (+ ipwadh_mata.do), ipwadh_eval, stcox_ipwadh (+ stcox_ipwadh_mata.do), ipwadh, ipwadh_mi (mi data)
* HOW TO RUN     cd to the folder that contains all .ado / .do files of the suite and the data/ folder (S3_N200, S8_N200), then
*                    do cscript_ipwadh.do                       (quick set: about the length of the individual tests)
*                    global ipwadh_cs_full 1                    (before the do: ALSO the slow N=2000 pipeline / scenario / multiple-imputation tests;
*                                                                these need data/S8_N2000 as well)
*                    global ipwadh_cs_only "test_ipwadh_eval"   (before the do: run one test only; blank = all)
* HOW A TEST COUNTS  A test passes when its do-file runs to the end (assert-style tests stop with an error at the first failure) AND, for tests that
*                    keep a counter, the global $nfail is 0.  Expected numbers inside the tests come from independent Python references (ref_*.py).
* NOT part of the certification: profile_ipwadh.do (timing only), the simulation scripts (*.py).
* mi data (Stata 16+): test_ipwadh_mi_wire (quick: ipwadh -> ipwadh_mi hand-over on 100 persons of S8_N200, M = 2) and, in the full set,
*                    test_ipwadh_mi_frame (ipwadh_mi.ado on 300 persons of S8_N2000, M = 3) and test_ipwadh_mi_strata (mi estimate: stcox with strata(), same data, M = 3).
*                    The non-mi path of ipwadh is certified by test_ipwadh_wrapper.
* ---------------------------------------------------------------------------------------------------------------------------------------
* -cscript- clears all global macros: keep the two switches in locals across it
local cs_full "$ipwadh_cs_full"
local cs_only `"$ipwadh_cs_only"'
cscript "ipwadh certification"
if "`cs_full'" != "" global ipwadh_cs_full "`cs_full'"
if `"`cs_only'"' != "" global ipwadh_cs_only `"`cs_only'"'
set more off
set linesize 150
capture log close _all
log using "cscript_ipwadh.log", replace text
if "$ipwadh_cs_full" == "" global ipwadh_cs_full 0
di as txt "ipwadh certification   Stata " c(stata_version) " " c(edition_real) "   " c(current_date) " " c(current_time)
di as txt "folder: `c(pwd)'"

* ----- 0. all files of the suite present? versions of the loaded commands ----------------------------------------------------------------------
local files "ipwadhc.ado ipwadh_covlists.ado ipwadh_split.ado ipwadh_mata.do ipwadh_eval.ado stcox_ipwadh.ado stcox_ipwadh_mata.do ipwadh.ado ipwadh_mi.ado"
local nmiss 0
di as txt _n "--- files and versions (the first line of each .ado: *! name version) ---"
foreach f of local files {
    capture findfile `f'
    if _rc {
        di as err "  MISSING  `f'"
        local ++nmiss
    }
    else {
        local fn `"`r(fn)'"'
        tempname fh
        file open `fh' using `"`fn'"', read text
        file read `fh' line
        file close `fh'
        di as txt "  " %-24s "`f'" as res substr(`"`line'"', 1, 70)
    }
}
foreach dd in data/S3_N200 data/S8_N200 {
    capture confirm file "`dd'/episodes.dta"
    if _rc {
        di as err "  MISSING  `dd'/episodes.dta"
        local ++nmiss
    }
}
if $ipwadh_cs_full == 1 {
    capture confirm file "data/S8_N2000/episodes.dta"
    if _rc {
        di as err "  MISSING  data/S8_N2000/episodes.dta  (needed by the full set)"
        local ++nmiss
    }
}
if `nmiss' > 0 {
    di as err _n "`nmiss' file(s) missing: the certification cannot run. Run it from the folder with all files and the data/ folder."
    log close
    exit 601
}

* ----- 1. the tests ------------------------------------------------------------------------------------------------------------------------------
local quick "test_ipwadh_split test_ipwadh_covlists test_ipwadhc test_ipwadhc_variants test_ipwadhc_terminal test_ipwadh_eval test_ipwadh_eval_trunc test_ipwadh_eval_warn test_stcox_ipwadh test_ipwadh_wrapper test_ipwadh_mi_wire"
local slow  "test_ipwadh_pipeline_N2000 test_ipwadh_scenarios_N2000 test_ipwadh_S8_mi test_ipwadh_mi_frame test_ipwadh_mi_strata"
local tests "`quick'"
if $ipwadh_cs_full == 1 local tests "`quick' `slow'"
if `"$ipwadh_cs_only"' != "" local tests `"$ipwadh_cs_only"'

local nt 0
local npass 0
local nfailed 0
foreach t of local tests {
    local ++nt
    local name`nt' "`t'"
    di as txt _n "{hline 150}" _n as res "RUNNING `t'.do" as txt "   (" c(current_time) ")" _n "{hline 150}"
    local t0 = clock(c(current_time), "hms")
    macro drop nfail nchk
    capture confirm file "`t'.do"
    if _rc {
        local rc`nt' = 601
        local ck`nt' "."
    }
    else {
        capture noisily do "`t'.do"
        local rc`nt' = _rc
        local ck`nt' "."
        if "$nchk" != "" local ck`nt' "$nchk"
        if `rc`nt'' == 0 & "$nfail" != "" & "$nfail" != "0" local rc`nt' = 9
    }
    local sec`nt' = mod(clock(c(current_time), "hms") - `t0' + 86400000, 86400000) / 1000
    if `rc`nt'' == 0 local ++npass
    else local ++nfailed
    * the test files start with -clear all-; make sure the certification itself keeps logging
    capture log query
    if _rc | "`r(status)'" != "on" {
        capture log using "cscript_ipwadh.log", append text
    }
}

* ----- 2. summary ---------------------------------------------------------------------------------------------------------------------------------
di as txt _n "{hline 100}" _n as res "ipwadh certification: summary" _n as txt "{hline 100}"
di as txt %-34s "test" %8s "result" %12s "counted" %10s "seconds"
forvalues i = 1/`nt' {
    local res = cond(`rc`i'' == 0, "pass", "FAIL r(" + string(`rc`i'') + ")")
    local col = cond(`rc`i'' == 0, "res", "err")
    di as txt %-34s "`name`i''" as `col' %8s "`res'" as txt %12s "`ck`i''" %10.0f `sec`i''
}
di as txt "{hline 100}"
di as txt "tests run: `nt'   passed: `npass'   failed: `nfailed'"
if `nfailed' == 0 {
    di as res "ALL TESTS PASSED"
}
else {
    di as err "CERTIFICATION FAILED: look at the FAIL lines in cscript_ipwadh.log"
}
if $ipwadh_cs_full != 1 & `"$ipwadh_cs_only"' == "" di as txt "(quick set only; global ipwadh_cs_full 1 adds the slow N=2000 pipeline, scenario and multiple-imputation tests)"
log close
if `nfailed' > 0 exit 9
