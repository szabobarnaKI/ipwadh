*! stcox_ipwadh 0.1.1  Weighted Cox / LWYY model with a robust variance that accounts for the ESTIMATION of the ipwadh weights  (NOT YET RUN IN STATA)
*! 0.1.1: compact default output (one header line, the stcox table, one line with the SE ratios); option detail restores the comparison table and the reproduction checks. Numbers unchanged. 0.1.0 passed test_stcox_ipwadh.do (57/59, 2 wrong expected counts since corrected).
*!
*! WHAT IT DOES
*!   Runs  stcox <varlist>, breslow vce(cluster id)  on data that were stset with the ipwadh weights, then replaces the cluster-robust
*!   variance (which treats the weights as known) by the variance of the stacked estimating equations: the weighted Cox score plus the
*!   scores of ALL Poisson rate models that produced the weights (adherence changes, death, loss to follow-up; numerator and denominator
*!   models, every stratum).  psi_i = U_i + sum_f B_f V_f S_if,  V = A^-1 (sum_i psi_i psi_i') A^-1 * G/(G-1)   (derivation: LOG_addendum_session5c.md,
*!   Python reference ref_se.py).  It is a SEPARATE command: ipwadhc / ipwadh_split / ipwadh_eval and a plain  stcox ..., vce(cluster id)  stay as they are.
*!
*! REQUIREMENTS (checked)
*!   1. ipwadhc (>= 0.3.1) was run with option  parts  on the unsplit episodes: it stores the model specification in the data (characteristics _dta[ipwadh_*]).
*!   2. ipwadh_split, then ipwadh_eval (creates _ipw_psw and, with truncate(), _ipw_pswt).  Data may NOT have been changed (rows, covariates) after ipwadhc.
*!   3. stset t1 [pw=_ipw_psw (or _ipw_pswt)] if adh==1, time0(t0) failure(...)    (pweights, no id(), no origin()/scale()).
*!   4. The Cox model contains arm plus exactly the stcoxvars() that were given to ipwadhc (strata() allowed).  Not allowed: efron/exactp ties, tvc(), shared(), offset.
*!   The rate models are REFITTED internally from the stored specification (cheap compared with the split); the recomputed weights are compared with the stored ones
*!   (r(maxdiff_parts), r(maxdiff_w)) and the recomputed cluster-robust variance with e(V) of stcox (r(maxdiff_V)): the command stops if they disagree.
*!
*! TRUNCATED WEIGHTS: pieces clipped by truncate() carry no derivative (the percentile bounds are treated as fixed); the variance of the bounds is ignored.
*!
*! SYNTAX
*!   stcox_ipwadh varlist [, STRata(varlist) NOHR Level(#) EVALat(end|mid|start) PREfix(name) GENerate(name) CLuster(varname) TOLerance(#) ]
*!   evalat() must be the value used in ipwadh_eval (default end).  prefix()/generate() as in ipwadh_eval (defaults: the prefix stored by ipwadhc).
*!   cluster(): default is the id() of ipwadhc.  tolerance(): largest allowed relative difference between recomputed and stored weights (default 1e-6).
*! RESULTS  e(V) = adjusted variance; e(V_naive) = the cluster-robust variance of stcox (weights treated as known); e(ipwadh_adjusted) = 1;
*!   e(ipwadh_maxdiff_parts), e(ipwadh_maxdiff_w), e(ipwadh_maxdiff_V), e(ipwadh_maxscore), e(ipwadh_nclipped).  All stcox results (e(b), e(N), ...) are unchanged.

program define stcox_ipwadh, eclass sortpreserve
    version 15
    syntax varlist(fv) [, STRata(varlist) NOHR Level(cilevel) EVALat(string) PREfix(name) GENerate(name) CLuster(varname) TOLerance(real 1e-6) DETail ]

    * (re)load the Mata functions every time (the file drops old versions first)
    quietly findfile stcox_ipwadh_mata.do
    quietly run `"`r(fn)'"'
    if "`evalat'" == "" local evalat end
    local evcode = cond("`evalat'" == "end", 1, cond("`evalat'" == "mid", 2, cond("`evalat'" == "start", 3, 0)))
    if `evcode' == 0 {
        di as err "evalat() must be end, mid or start"
        exit 198
    }

    * ---------- stset data with pweights ----------
    capture quietly st_is 2 analysis
    if _rc {
        di as err "data not st"
        exit 119
    }
    local wt : char _dta[st_wt]
    local wv : char _dta[st_wv]
    local stid : char _dta[st_id]
    if "`wt'" != "pweight" {
        di as err "the data must be stset with pweights:  stset t1 [pw=_ipw_psw (or _ipw_pswt)] if adh==1, time0(t0) failure(...)"
        exit 198
    }
    if "`stid'" != "" {
        di as err "stset with id() is not supported (Stata refuses time-varying weights with id()); stset without id() and without exit()"
        exit 198
    }

    * ---------- specification stored by ipwadhc ----------
    foreach c in id time0 time1 adh arm cens byarm rn_a rd_a rn_c rd_c popts prefix {
        local sp_`c' : char _dta[ipwadh_`c']
    }
    if "`sp_id'" == "" {
        di as err "no ipwadhc specification found in the data (characteristics _dta[ipwadh_*])."
        di as err "Run ipwadhc (version 0.3.1 or later) with option parts on the unsplit episodes, then ipwadh_split and ipwadh_eval, and do not drop or change rows afterwards."
        exit 459
    }
    if "`prefix'" == "" local prefix `sp_prefix'
    if "`generate'" == "" local generate `prefix'
    if "`cluster'" == "" local cluster `sp_id'
    foreach v in `prefix'_orig `prefix'_piece `prefix'_E0 `prefix'_J0 `prefix'_dA `prefix'_Ec0 `prefix'_dC `generate'_psw `sp_time0' `sp_time1' `sp_id' `sp_adh' `sp_arm' `cluster' {
        capture confirm numeric variable `v'
        if _rc {
            di as err "variable `v' not found (or not numeric): run ipwadhc ... parts, ipwadh_split, ipwadh_eval first"
            exit 111
        }
    }
    foreach v of local sp_cens {
        confirm numeric variable `v'
    }
    if "`wv'" != "`generate'_psw" & "`wv'" != "`generate'_pswt" {
        di as err "the stset weight is `wv'; it must be `generate'_psw (untruncated) or `generate'_pswt (truncated) made by ipwadh_eval"
        exit 198
    }
    local ncens : word count `sp_cens'

    * ---------- the Cox model ----------
    local runopts breslow vce(cluster `cluster')
    if "`strata'" != "" local runopts `runopts' strata(`strata')
    capture estimates drop __ipwse_cox
    quietly stcox `varlist', `runopts'
    quietly estimates store __ipwse_cox
    if "`e(cmd)'" != "cox" & "`e(cmd)'" != "stcox" {
        di as err "stcox did not run as expected (e(cmd) = `e(cmd)')"
        exit 198
    }
    * unsupported extras (tvc/texp, offset, shared) are caught by the built-in check that the recomputed cluster-robust variance equals e(V)
    tempvar esamp
    quietly gen byte `esamp' = e(sample)
    tempname Vcox Vadj Vnaive
    matrix __ipwse_b = e(b)
    matrix `Vcox' = e(V)
    matrix __ipwse_V = e(V)
    _ms_omit_info __ipwse_b
    matrix __ipwse_om = r(omit)
    local xl ""
    local cn : colnames __ipwse_b
    local k 0
    foreach c of local cn {
        local ++k
        if __ipwse_om[1, `k'] == 1 continue
        fvrevar `c'
        local xl `xl' `r(varlist)'
    }
    local es "`e(strata)'"
    tempvar cstr
    if "`es'" != "" quietly egen long `cstr' = group(`es') if `esamp'
    else quietly gen byte `cstr' = 1

    * ---------- row-level variables (full data set) ----------
    local t0 `sp_time0'
    local t1 `sp_time1'
    tempvar rlen rowsel chg sa sc
    quietly bysort `prefix'_orig (`prefix'_piece): gen double `rlen' = `t1'[_N] - `t0'[1]
    local ycs ""
    forvalues c = 1/`ncens' {
        local cv : word `c' of `sp_cens'
        tempvar yc`c'
        quietly bysort `prefix'_orig (`prefix'_piece): gen byte `yc`c'' = `cv'[_N]
        local ycs `ycs' `yc`c''
    }
    quietly gen byte `rowsel' = (`prefix'_piece == 1) & !missing(`prefix'_E0)

    * ---------- refit the rate models on the original rows (reduced data set) ----------
    preserve
    quietly keep if `rowsel'
    quietly sort `sp_id' `t0'
    quietly by `sp_id': gen byte `chg' = cond(_n < _N, `sp_adh' != `sp_adh'[_n+1], 0)
    if "`sp_byarm'" == "1" {
        quietly egen long `sa' = group(`sp_arm' `sp_adh')
        quietly egen long `sc' = group(`sp_arm')
    }
    else {
        quietly gen double `sa' = `sp_adh'
        quietly gen double `sc' = 1
    }
    mata: ipwadh_se_start()
    mata: ipwadh_se_rows("`sp_id'", "`t0'", "`rlen'", "`sp_adh'", "`prefix'_orig", "`chg'", "`prefix'_E0 `prefix'_J0 `prefix'_dA `prefix'_Ec0 `prefix'_dC", `ncens')
    local f 0
    quietly levelsof `sa', local(lva)
    foreach l of local lva {
        local ++f
        _stcox_ipwadh_fit `f' 1 1 0 `chg' `"`sp_rn_a'"' `sa' `l' `rlen' `"`sp_popts'"' "adherence changes, reduced model, stratum `l'"
        local ++f
        _stcox_ipwadh_fit `f' 1 2 0 `chg' `"`sp_rd_a'"' `sa' `l' `rlen' `"`sp_popts'"' "adherence changes, full model, stratum `l'"
    }
    quietly levelsof `sc', local(lvc)
    forvalues c = 1/`ncens' {
        local yv : word `c' of `ycs'
        local cl : word `c' of `sp_cens'
        foreach l of local lvc {
            local ++f
            _stcox_ipwadh_fit `f' 2 1 `c' `yv' `"`sp_rn_c'"' `sc' `l' `rlen' `"`sp_popts'"' "`cl', reduced model, stratum `l'"
            local ++f
            _stcox_ipwadh_fit `f' 2 2 `c' `yv' `"`sp_rd_c'"' `sc' `l' `rlen' `"`sp_popts'"' "`cl', full model, stratum `l'"
        }
    }
    restore
    * the Poisson refits overwrote e(): bring back the stcox results
    quietly estimates restore __ipwse_cox

    * ---------- the variance ----------
    capture noisily mata: ipwadh_se_final("_t0", "_t", "_d", "`wv'", "`generate'_psw", "`prefix'_orig", "`cstr'", "`esamp'", "`t0'", "`t1'", "`xl'", "__ipwse_b", "__ipwse_om", "__ipwse_V", "__ipwse_Vadj", `evcode', 1)
    local rc = _rc
    if !`rc' {
        * st_matrix() in Mata drops the row/column names: restore them from e(V)
        local fn : colfullnames `Vcox'
        capture matrix colnames __ipwse_Vadj = `fn'
        local rc = _rc
        if !`rc' matrix rownames __ipwse_Vadj = `fn'
        if `rc' di as err "stcox_ipwadh: dimension mismatch between the adjusted variance (" rowsof(__ipwse_Vadj) "x" colsof(__ipwse_Vadj) ") and e(V) (" rowsof(`Vcox') "x" colsof(`Vcox') ")"
    }
    capture mata: ipwadh_se_clear()
    if `rc' {
        di as err "stcox_ipwadh: computation failed (r(`rc'))"
        exit `rc'
    }
    local maxdP = scalar(__ipwse_maxdP)
    local maxdw = scalar(__ipwse_maxdw)
    local maxdV = scalar(__ipwse_maxdV)
    local maxU = scalar(__ipwse_maxU)
    local Gc = scalar(__ipwse_G)
    local nclip = scalar(__ipwse_nclip)
    if `maxdP' > `tolerance' | `maxdw' > `tolerance' {
        di as err "stcox_ipwadh: the refitted rate models do not reproduce the stored weights (max relative difference: parts " %9.2e `maxdP' ", weights " %9.2e `maxdw' ")."
        di as err "  Possible reasons: rows or covariates changed after ipwadhc; evalat() differs from the ipwadh_eval call; different Poisson options; tolerance() too small."
        exit 498
    }
    if `maxdV' > 1e-3 {
        di as err "stcox_ipwadh: the recomputed cluster-robust variance differs from the e(V) of stcox (max relative difference of the diagonal " %9.2e `maxdV' ")."
        di as err "  The weighted Cox model is not the one assumed (ties other than Breslow? stset origin()/scale()? weights?). The adjusted variance is NOT reported."
        exit 498
    }

    * ---------- post ----------
    ereturn repost V = __ipwse_Vadj, esample(`esamp')
    matrix `Vnaive' = `Vcox'
    ereturn matrix V_naive = `Vnaive'
    ereturn scalar ipwadh_adjusted = 1
    ereturn scalar ipwadh_maxdiff_parts = `maxdP'
    ereturn scalar ipwadh_maxdiff_w = `maxdw'
    ereturn scalar ipwadh_maxdiff_V = `maxdV'
    ereturn scalar ipwadh_maxscore = `maxU'
    ereturn scalar ipwadh_nclipped = `nclip'
    capture estimates drop __ipwse_cox
    ereturn local ipwadh_vce "stacked estimating equations: Cox score + Poisson weight models"

    if "`level'" == "" local level 95
    di
    di as txt "stcox_ipwadh: standard errors adjusted for the estimated weights (" as res `f' as txt " Poisson fits refitted); cluster-robust variance in e(V_naive)"
    stcox, `nohr' level(`level')
    if "`detail'" != "" {
        di
        di as txt "Standard errors of the coefficients" _col(48) "cluster-robust" _col(66) "adjusted" _col(80) "ratio"
        di as txt "{hline 86}"
    }
    local ratios ""
    local k 0
    foreach c of local cn {
        local ++k
        if __ipwse_om[1, `k'] == 1 continue
        local sn = sqrt(el(e(V_naive), `k', `k'))
        local sa2 = sqrt(el(e(V), `k', `k'))
        if "`detail'" != "" di as txt "  " abbrev("`c'", 40) _col(48) as res %12.6f `sn' _col(62) %12.6f `sa2' _col(78) %7.4f `sa2' / `sn'
        local rr = string(`sa2' / `sn', "%6.3f")
        local ratios `"`ratios'   `c' `rr'"'
    }
    if "`detail'" != "" {
        di as txt "{hline 86}"
        di as txt "  cluster-robust = weights treated as known (plain stcox, vce(cluster)); adjusted = also accounts for the estimated rate models."
    }
    else di as txt "SE ratio adjusted / cluster-robust:" as res `"`ratios'"'
    if `nclip' > 0 di as txt "  note: " as res `nclip' as txt " pieces carry truncated weights (bounds treated as fixed: no weight derivative)."
    if "`detail'" != "" di as txt "  checks: weights reproduced to " as res %9.2e max(`maxdP', `maxdw') as txt ", cluster-robust variance reproduced to " as res %9.2e `maxdV' as txt ", max |score| at beta-hat " as res %9.2e `maxU'

    foreach m in __ipwse_b __ipwse_V __ipwse_om __ipwse_Vadj __ipwse_Vadj_naive __ipwse_sb __ipwse_sV __ipwse_so {
        capture matrix drop `m'
    }
    foreach sc in __ipwse_maxdP __ipwse_maxdw __ipwse_maxdV __ipwse_maxU __ipwse_G __ipwse_nclip __ipwse_npieces {
        capture scalar drop `sc'
    }
end


program define _stcox_ipwadh_fit
    args f fam role cause dep rhs strat lvl len popts what
    tempvar m
    quietly gen byte `m' = (`strat' == `lvl')
    quietly count if `m' & `dep' == 1
    if r(N) == 0 {
        di as err "no events for the model of `what'"
        exit 2000
    }
    local po `"`popts'"'
    local usedtight 0
    if `"`po'"' == "" {
        local po "tolerance(1e-12) ltolerance(1e-14) nrtolerance(1e-12)"
        local usedtight 1
    }
    capture quietly poisson `dep' `rhs' if `m', exposure(`len') `po'
    if `usedtight' & (_rc | (_rc == 0 & e(converged) == 0)) {
        capture quietly poisson `dep' `rhs' if `m', exposure(`len')
    }
    if _rc {
        local rc = _rc
        di as err "Poisson model for `what' failed (r(`rc'))"
        exit `rc'
    }
    matrix __ipwse_sb = e(b)
    matrix __ipwse_sV = e(V)
    _ms_omit_info __ipwse_sb
    matrix __ipwse_so = r(omit)
    local zl ""
    local cn : colnames __ipwse_sb
    local k 0
    foreach c of local cn {
        local ++k
        if "`c'" == "_cons" continue
        if __ipwse_so[1, `k'] == 1 continue
        fvrevar `c'
        local zl `zl' `r(varlist)'
    }
    mata: ipwadh_se_addfit(`f', `fam', `role', `cause', "`zl'", "`m'", "`dep'", "__ipwse_sb", "__ipwse_so", "__ipwse_sV")
end
