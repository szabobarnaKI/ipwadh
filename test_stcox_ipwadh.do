* test_stcox_ipwadh.do -- run in Stata 19 from the folder containing ipwadhc.ado (v0.3.1), ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do,
* ipwadh_eval.ado (v0.3.0), stcox_ipwadh.ado, stcox_ipwadh_mata.do and data/.   NOT YET RUN.   Data: S3_N200 (seed 1), exact split (77,718 pieces), evalat(end).
* Expected values: independent Python reference ref_se.py (se_expected.py -> se_expected_S3_N200.json).  Its influence function was verified against a
* full-pipeline numerical perturbation (test_if_ref.py) and its derivative terms against finite differences (test_se_ref.py).
* Cox model: i.arm c.age i.male (arm + exactly the stcoxvars given to ipwadhc).  "naive" = cluster-robust SE of plain stcox (weights treated as known).
clear all
set more off
global nfail = 0
global nchk = 0
local d "data/S3_N200"
local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm)"
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
capture program drop chkr
program define chkr
    args label value expected tol
    global nchk = $nchk + 1
    if abs(`value' / `expected' - 1) < `tol' di as txt "  PASS  `label': " %14.8f `value' "  (expected " %14.8f `expected' ", rel. tol `tol')"
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label': " %14.8f `value' "  (expected " %14.8f `expected' ", rel. tol `tol')"
    }
end
capture program drop sechk
program define sechk
    args cname eb en ea
    tempname Vn Va
    matrix `Vn' = e(V_naive)
    matrix `Va' = e(V)
    local j = colnumb(`Va', "`cname'")
    local b = _b[`cname']
    local sn = sqrt(`Vn'[`j', `j'])
    local sa = sqrt(`Va'[`j', `j'])
    chk  "beta `cname'" `b' `eb' 1e-5
    chkr "SE cluster-robust (naive) `cname'" `sn' `en' 2e-4
    chkr "SE adjusted `cname'" `sa' `ea' 2e-4
end
capture program drop expect_rc
program define expect_rc
    args label rc expected
    global nchk = $nchk + 1
    if `rc' == `expected' di as txt "  PASS  `label': r(`rc')"
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label': r(`rc') (expected r(`expected'))"
    }
end

* ===== A. death model only ==========================================================================================================================
use "`d'/episodes.dta", clear
ipwadhc, `base' dead(dead) `lists' parts replace
ipwadh_split, `sp'
assert _N == 77718
ipwadh_eval, time0(t0) time1(t1) adherence(adh) truncate(1 99) replace

di as res _n "=== 1. death only, untruncated weights ==="
quietly stset t1 [pw=_ipw_psw] if adh==1, time0(t0) failure(inf)
quietly stcox i.arm c.age i.male, breslow vce(cluster id)
matrix Vplain = e(V)
stcox_ipwadh i.arm c.age i.male
sechk 1.arm  -0.3065022  0.15855301 0.15800579
sechk age     0.03085327 0.00564575 0.00567065
sechk 1.male  0.31727187 0.14849424 0.14830639
chk "G (clusters)" e(N_clust) 192 0.5
chk "no truncated pieces" e(ipwadh_nclipped) 0 0.5
local ok1 = (e(ipwadh_maxdiff_parts) < 1e-6)
local ok2 = (e(ipwadh_maxdiff_w) < 1e-6)
local ok3 = (e(ipwadh_maxdiff_V) < 1e-4)
chk "refit reproduces stored parts (<1e-6)" `ok1' 1 0.5
chk "refit reproduces weights (<1e-6)" `ok2' 1 0.5
chk "recomputed cluster-robust V = stcox V (<1e-4)" `ok3' 1 0.5
di as txt "  info: max rel. diff parts " %9.2e e(ipwadh_maxdiff_parts) ", weights " %9.2e e(ipwadh_maxdiff_w) ", V " %9.2e e(ipwadh_maxdiff_V) ", max |score| " %9.2e e(ipwadh_maxscore)
* e(V_naive) is the plain stcox variance
tempname Vn
matrix `Vn' = e(V_naive)
local dmax = 0
forvalues i = 1/`=rowsof(Vplain)' {
    local dd = abs(`Vn'[`i', `i'] - Vplain[`i', `i'])
    if `dd' > `dmax' local dmax = `dd'
}
chk "e(V_naive) == plain stcox e(V) (max abs diff of the diagonal)" `dmax' 0 1e-12
* estimates are those of stcox
quietly stcox i.arm c.age i.male, breslow vce(cluster id)
local b_plain = _b[1.arm]
quietly stcox_ipwadh i.arm c.age i.male
chk "e(b) unchanged" _b[1.arm] `b_plain' 1e-12

di as res _n "=== 2. death only, truncated weights p1/p99 (pieces at the bounds have no weight derivative) ==="
quietly stset t1 [pw=_ipw_pswt] if adh==1, time0(t0) failure(inf)
stcox_ipwadh i.arm c.age i.male
sechk 1.arm  -0.29944968 0.16017776 0.16256843
sechk age     0.0296949  0.00550962 0.00567031
sechk 1.male  0.33469784 0.1483709  0.14967054
* (Python counted 1384 clipped pieces over ALL pieces; Stata counts the pieces of the Cox sample only, i.e. adherent ones: 1294)
chk "clipped pieces in the Cox sample" e(ipwadh_nclipped) 1294 0.5

di as res _n "=== 3. stratified Cox (strata(male)), untruncated ==="
quietly stset t1 [pw=_ipw_psw] if adh==1, time0(t0) failure(inf)
stcox_ipwadh i.arm c.age, strata(male)
sechk 1.arm  -0.28346841 0.15699487 0.15636029
sechk age     0.03078317 0.00591841 0.00593319

di as res _n "=== 4. option checks: nohr, level(90), evalat(mid) must fail the weight reproduction ==="
stcox_ipwadh i.arm c.age i.male, nohr level(90)
capture noisily stcox_ipwadh i.arm c.age i.male, evalat(mid)
expect_rc "evalat(mid) does not reproduce the stored end-weights" _rc 498
capture noisily stcox_ipwadh i.arm c.age i.male, evalat(sideways)
expect_rc "bad evalat()" _rc 198

di as res _n "=== 5. error cases ==="
quietly stset t1 [iw=_ipw_psw] if adh==1, time0(t0) failure(inf)
capture noisily stcox_ipwadh i.arm c.age i.male
expect_rc "iweights" _rc 198
quietly stset t1 [pw=_ipw_sw] if adh==1, time0(t0) failure(inf)
capture noisily stcox_ipwadh i.arm c.age i.male
expect_rc "weight variable is the row-end weight, not _ipw_psw/_ipw_pswt" _rc 198
quietly stset t1 [pw=_ipw_psw] if adh==1, time0(t0) failure(inf)
preserve
replace age = age + 3 * (mod(id, 7) == 0)
capture noisily stcox_ipwadh i.arm c.age i.male
expect_rc "covariate changed after ipwadhc" _rc 498
restore
preserve
char _dta[ipwadh_id]
capture noisily stcox_ipwadh i.arm c.age i.male
expect_rc "no specification stored" _rc 459
restore
preserve
capture drop _ipw_E0
capture noisily stcox_ipwadh i.arm c.age i.male
expect_rc "parts missing" _rc 111
restore
capture noisily stcox_ipwadh i.arm c.age i.male, efron
expect_rc "efron is not an option of stcox_ipwadh" _rc 198

* ===== B. death and loss to follow-up ================================================================================================================
di as res _n "=== 6. death + ltfu models, untruncated and truncated ==="
use "`d'/episodes.dta", clear
ipwadhc, `base' dead(dead) ltfu(ltfu) `lists' parts replace
ipwadh_split, `sp'
ipwadh_eval, time0(t0) time1(t1) adherence(adh) truncate(1 99) replace
quietly stset t1 [pw=_ipw_psw] if adh==1, time0(t0) failure(inf)
stcox_ipwadh i.arm c.age i.male
sechk 1.arm  -0.31047101 0.15863278 0.15800211
sechk age     0.03074292 0.00560891 0.00563517
sechk 1.male  0.31446247 0.14863171 0.14855744
quietly stset t1 [pw=_ipw_pswt] if adh==1, time0(t0) failure(inf)
stcox_ipwadh i.arm c.age i.male
sechk 1.arm  -0.30399958 0.16048686 0.16261854
sechk age     0.02952933 0.00548937 0.00564359
sechk 1.male  0.33346846 0.14868973 0.14996666
chk "clipped pieces in the Cox sample" e(ipwadh_nclipped) 1286 0.5

di as res _n "=== summary ==="
di as txt "checks: " $nchk "   failed: " $nfail
if $nfail == 0 di as res "ALL PASSED"
else di as err "SOME FAILED"
