* test_ipwadh_wrapper.do -- run in Stata 19 from the folder with ipwadhc.ado (0.3.1), ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do,
* ipwadh_eval.ado, stcox_ipwadh.ado, stcox_ipwadh_mata.do, ipwadh.ado (0.1.1) and data/.   NOT YET RUN.   Data: S3_N200, exact split.
* Expected values: the same Python reference values as test_stcox_ipwadh.do (death model only).
clear all
set more off
global nfail = 0
global nchk = 0
local d "data/S3_N200"
local common "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) failure(inf) dead(dead) stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"

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

di as res _n "=== 1. default: truncate(1 99), cluster-robust SE; original data restored ==="
use "`d'/episodes.dta", clear
local N0 = _N
quietly ds
local vars0 "`r(varlist)'"
ipwadh, `common'
chk "beta 1.arm (truncated)" _b[1.arm] -0.29944968 1e-5
chk "SE cluster-robust 1.arm" _se[1.arm] 0.16017776 2e-5
chk "beta age" _b[age] 0.0296949 1e-5
chk "beta 1.male" _b[1.male] 0.33469784 1e-5
chk "G clusters" e(N_clust) 192 0.5
chk "data restored: N" _N `N0' 0.5
quietly ds
local same = ("`r(varlist)'" == "`vars0'")
chk "data restored: variables" `same' 1 0.5
local t1 = ("`e(ipwadh_se)'" == "cluster-robust")
chk "e(ipwadh_se) is cluster-robust" `t1' 1 0.5
local t2 = ("`e(ipwadh_weight)'" == "_ipw_pswt")
chk "e(ipwadh_weight) is _ipw_pswt" `t2' 1 0.5
chk "e(ipwadh_n_pieces)" e(ipwadh_n_pieces) 77718 0.5
chk "ESS truncated about 51%" e(ipwadh_ess_trunc) 0.510 0.005
chk "ESS untruncated about 45.6%" e(ipwadh_ess) 0.456 0.005

di as res _n "=== 2. notruncate ==="
use "`d'/episodes.dta", clear
ipwadh, `common' notruncate
chk "beta 1.arm (untruncated)" _b[1.arm] -0.30650220 1e-5
chk "SE cluster-robust 1.arm" _se[1.arm] 0.15855301 2e-5
local t3 = ("`e(ipwadh_weight)'" == "_ipw_psw")
chk "e(ipwadh_weight) is _ipw_psw" `t3' 1 0.5

di as res _n "=== 3. adjustse (truncated and untruncated) ==="
use "`d'/episodes.dta", clear
ipwadh, `common' adjustse
chk "beta 1.arm" _b[1.arm] -0.29944968 1e-5
chk "SE adjusted 1.arm (truncated)" _se[1.arm] 0.16256843 3e-5
matrix Vn = e(V_naive)
local sn = sqrt(Vn[2,2])
chk "SE naive from e(V_naive)" `sn' 0.16017776 3e-5
local t4 = ("`e(ipwadh_se)'" == "adjusted")
chk "e(ipwadh_se) is adjusted" `t4' 1 0.5
use "`d'/episodes.dta", clear
ipwadh, `common' adjustse notruncate
chk "SE adjusted 1.arm (untruncated)" _se[1.arm] 0.15800579 3e-5

di as res _n "=== 4. keepdata: split data, weights and stset stay; same result as the manual pipeline ==="
use "`d'/episodes.dta", clear
ipwadh, `common' keepdata
local b1 = _b[1.arm]
chk "split data kept" _N 77718 0.5
confirm variable _ipw_pswt
confirm variable _ipw_psw
quietly stcox i.arm c.age i.male, breslow vce(cluster id)
chk "manual stcox on the kept data gives the same beta" _b[1.arm] `b1' 1e-12
capture noisily ipwadh, `common'
expect_rc "second call on split data is refused" `=_rc' 110

di as res _n "=== 5. weightsonly (implies keepdata), then own stcox ==="
use "`d'/episodes.dta", clear
ipwadh, `common' weightsonly
chk "split data kept (weightsonly)" _N 77718 0.5
quietly stcox i.arm c.age i.male, breslow vce(cluster id)
chk "beta 1.arm" _b[1.arm] -0.29944968 1e-5

di as res _n "=== 6. strata(male): male removed from the covariates, kept in the weights ==="
use "`d'/episodes.dta", clear
ipwadh, `common' strata(male) notruncate
chk "beta 1.arm (stratified, untruncated)" _b[1.arm] -0.28346841 1e-5
chk "SE 1.arm (stratified)" _se[1.arm] 0.15699487 3e-5
local t5 = ("`e(strata)'" != "")
chk "stratified model (e(strata) set)" `t5' 1 0.5

di as res _n "=== 7. if/in: only the selected persons are used; the rest of the data is unchanged ==="
use "`d'/episodes.dta", clear
local N0 = _N
ipwadh if id <= 100, `common'
chk "data restored (if)" _N `N0' 0.5
chk "only persons with id <= 100 (96 have episodes)" e(N_clust) 96 0.5

di as res _n "=== 8. option errors ==="
use "`d'/episodes.dta", clear
capture ipwadh, `common' truncate(1 99) notruncate
expect_rc "truncate + notruncate" `=_rc' 198
capture ipwadh, `common' weightsonly adjustse
expect_rc "weightsonly + adjustse" `=_rc' 198
capture ipwadh, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) failure(inf) dead(dead) stcoxvars(c.age i.male) truncate(1)
local rcx = (_rc != 0)
expect_rc "truncate with one number is refused" `rcx' 1

di as res _n "=== 9. verbose (detail output of every step) gives the same result; detail options of the single commands run ==="
use "`d'/episodes.dta", clear
ipwadh, `common' verbose adjustse
chk "beta 1.arm (verbose, adjustse)" _b[1.arm] -0.29944968 1e-5
chk "SE adjusted (verbose)" _se[1.arm] 0.16256843 3e-5
use "`d'/episodes.dta", clear
ipwadhc, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker) parts replace detail
ipwadh_split, id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu) detail
ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) truncate(1 99) replace detail
quietly stset t1 [pw=_ipw_pswt] if adh==1, time0(t0) failure(inf)
stcox_ipwadh i.arm c.age i.male, detail
chk "beta 1.arm (single commands, detail)" _b[1.arm] -0.29944968 1e-5

di as res _n "=== summary ==="
di as txt "checks: $nchk   failed: $nfail"
if $nfail == 0 di as res "ALL PASSED"
else di as err "SOME CHECKS FAILED"
