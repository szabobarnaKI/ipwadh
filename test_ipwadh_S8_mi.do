* test_ipwadh_S8_mi.do -- scenario S8: MAR missingness in biomarker and abx (about 15% of rows each, 95% of persons have at least one missing row)
* + multiple imputation (mi impute chained, M = 10) + the ipwadh pipeline run INSIDE each imputation ("within" approach) + Rubin's rules.  NOT YET RUN.
* Run in Stata 19 from the folder containing ipwadhc.ado (v0.3.0), ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do, ipwadh_eval.ado (v0.2.0)
* and data/S8_N2000/episodes.dta (send file from this session; or regenerate: python3 ipwadh_simdata.py --scenario S8 --n 2000 --seed 20261005 --out data).
* RUNTIME: one pass of ipwadhc + exact split (7.5M rows) + 2 x stcox is about 4-6 min -> 40-60 min for M = 10. Set  local M 3  for a trial, or  local bandopt "bandpoints(auto)"
*          (bands: approximate, see ipwadh_split) to make it faster.  RAM about 2 GB.
* WHY NOT  mi estimate: ... : the pipeline SPLITS the data (adds rows) inside every imputation; mi estimate needs a command that leaves the data unchanged.
*   So the loop below does  mi extract m  -> pipeline -> store b and its robust variance -> Rubin's rules (qbar, W, B, T = W + (1+1/M) B, df, FMI) by hand.
* BENCHMARK = complete data: S8 is EXACTLY the S3 data set with values set to missing (verified in Python: all other columns identical), so the complete-data results of
*   the same pipeline on S3_N2000 (death + ltfu IPCW, exact split, evalat end; Python reference, Stata cluster SE) are the target of the MI analysis:
*     untruncated  coef -0.391823  SE 0.053725        truncated p1/p99  coef -0.390331  SE 0.050138
*   (marginal truth about -0.43). The imputation is random, so there is no digit-by-digit reference: the checks below are soft (PASS/FAIL, tolerance on the difference to the benchmark).
* Imputation model: biomarker (regress) and abx (logit, augment) on arm age male adh i.tsi ninf inf dead ltfu na t1; na = Nelson-Aalen cumulative hazard of infection (pooled over rows,
*   evaluated at the row end: a stand-in for the cumulative-hazard term recommended for imputation with a time-to-event outcome). Rows, not persons, are imputed (the data are
*   episodes); the person-level correlation of the biomarker (carried forward between visits) is NOT modelled: a limitation to keep in mind (the weights only need it as a covariate).
clear all
set more off
global nfail = 0
global nchk = 0
local M 10
local bandopt ""
local seed 20261010

capture program drop chk
program define chk
    args label value expected tol mode
    if "`mode'" == "" local mode "abs"
    if "`mode'" == "abs"  local dif = abs(`value' - `expected')
    if "`mode'" == "rel2" local dif = abs(`value' - `expected') / abs(`expected')
    global nchk = $nchk + 1
    if `dif' < `tol' {
        di as txt "  PASS  `label': " %12.6f `value' "  (benchmark " %12.6f `expected' ")"
    }
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label': " %12.6f `value' "  (benchmark " %12.6f `expected' ", `mode' diff " %10.3e `dif' ")"
    }
end

local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) ltfu(ltfu)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
timer clear
timer on 1

* ===== 0. the data and the documented behaviour on incomplete data ===================================================================
use "data/S8_N2000/episodes.dta", clear
assert _N == 69645
quietly count if missing(biomarker)
di as txt "rows with missing biomarker: " r(N) "   missing abx: " _continue
quietly count if missing(abx)
di r(N)
capture noisily ipwadhc, `base' `lists' replace
di as txt "ipwadhc on the incomplete data: rc = " _rc "  (expected 416: missing values in the weight-model covariates are reported, nothing is silently dropped)"
global nchk = $nchk + 1
if _rc != 416 {
    global nfail = $nfail + 1
    di as err "  FAIL  expected rc 416 on incomplete data, got " _rc
}
else di as txt "  PASS  rc 416 on incomplete data"

* ===== 1. imputation ==================================================================================================================
stset t1, time0(t0) failure(inf)
quietly sts generate na = na
stset, clear
quietly replace na = 0 if missing(na)      // zero-length rows (ignored by stset) have no Nelson-Aalen value
mi set wide
mi register imputed biomarker abx
mi register regular arm age male adh tsi ninf inf dead ltfu na t0 t1
mi impute chained (regress) biomarker (logit, augment) abx = arm age male adh i.tsi ninf inf dead ltfu na t1, add(`M') rseed(`seed')
save "S8_mi_tmp.dta", replace
di as res "1. imputation done: M = `M'"

* ===== 2. pipeline inside every imputation =================================================================================================
matrix bm = J(`M', 1, .)
matrix vm = J(`M', 1, .)
matrix bt = J(`M', 1, .)
matrix vt = J(`M', 1, .)
forvalues m = 1/`M' {
    di as res _n "--- imputation `m' of `M' ---"
    use "S8_mi_tmp.dta", clear
    mi extract `m', clear
    quietly count if missing(biomarker) | missing(abx)
    assert r(N) == 0
    quietly ipwadhc, `base' `lists' parts replace
    quietly summarize _ipw_sw if t1 > t0
    di as txt "  row-end weights: mean " %7.4f r(mean) "  max " %8.3f r(max)
    quietly ipwadh_split, id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu) `bandopt'
    quietly ipwadh_eval, time0(t0) time1(t1) adherence(adh) truncate(1 99) replace
    quietly stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
    quietly stcox i.arm, breslow vce(cluster id)
    matrix bm[`m', 1] = _b[1.arm]
    matrix vm[`m', 1] = _se[1.arm]^2
    quietly stset t1 [pw=_ipw_pswt] if adh == 1, time0(t0) failure(inf)
    quietly stcox i.arm, breslow vce(cluster id)
    matrix bt[`m', 1] = _b[1.arm]
    matrix vt[`m', 1] = _se[1.arm]^2
    di as txt "  untruncated coef " %9.5f bm[`m', 1] "  SE " %8.5f sqrt(vm[`m', 1]) "     truncated coef " %9.5f bt[`m', 1] "  SE " %8.5f sqrt(vt[`m', 1])
}
erase "S8_mi_tmp.dta"

* ===== 3. Rubin's rules ============================================================================================================================
capture mata: mata drop ipwadh_rubin()
mata:
void ipwadh_rubin(string scalar bn, string scalar vn, string scalar tag)
{
    real colvector b, v
    real scalar M, qbar, W, B, T, df, r, fmi
    b = st_matrix(bn); v = st_matrix(vn)
    M = rows(b)
    qbar = mean(b); W = mean(v); B = variance(b)
    T = W + (1 + 1/M) * B
    r = (1 + 1/M) * B / W
    df = (B > 0) ? (M - 1) * (1 + 1/r)^2 : .
    fmi = (B > 0) ? (r + 2/(df + 3)) / (r + 1) : 0
    st_numscalar(tag + "_q", qbar)
    st_numscalar(tag + "_se", sqrt(T))
    st_numscalar(tag + "_df", df)
    st_numscalar(tag + "_fmi", fmi)
    st_numscalar(tag + "_W", W)
    st_numscalar(tag + "_B", B)
}
end
mata: ipwadh_rubin("bm", "vm", "ru")
mata: ipwadh_rubin("bt", "vt", "rt")
foreach t in ru rt {
    local lab = cond("`t'" == "ru", "untruncated", "truncated p1/p99")
    local tq = cond(missing(scalar(`t'_df)), invnormal(0.975), invttail(scalar(`t'_df), 0.025))
    di as res _n "Rubin's rules, `lab' (M = `M'):"
    di as txt "  log HR " %9.5f scalar(`t'_q) "   SE " %8.5f scalar(`t'_se) "   95% CI " %8.5f (scalar(`t'_q) - `tq' * scalar(`t'_se)) " to " %8.5f (scalar(`t'_q) + `tq' * scalar(`t'_se))
    di as txt "  HR " %7.4f exp(scalar(`t'_q)) "   within-imputation variance " %10.3e scalar(`t'_W) "   between " %10.3e scalar(`t'_B) "   df " %8.1f scalar(`t'_df) "   FMI " %6.3f scalar(`t'_fmi)
}

* ===== 4. soft checks against the complete-data benchmark (S3 = S8 before values were set to missing) ============================================
di as res _n "=== 4. checks against the complete-data benchmark ==="
chk "MI untruncated log HR" scalar(ru_q) -0.391823 0.04
chk "MI untruncated SE (Rubin)" scalar(ru_se) 0.053725 0.20 rel2
chk "MI truncated log HR" scalar(rt_q) -0.390331 0.04
chk "MI truncated SE (Rubin)" scalar(rt_se) 0.050138 0.20 rel2
local nok = 0
forvalues m = 1/`M' {
    if !missing(bm[`m', 1]) local ++nok
}
chk "imputations with a result" `nok' `M' 0.5

timer off 1
timer list 1
di as res _n "test_ipwadh_S8_mi.do: $nchk checks, $nfail failed"
di as txt "Send me the whole log (the per-imputation lines and the Rubin tables are the result)."
