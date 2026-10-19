* test_ipwadh_scenarios_N2000.do -- run in Stata 19 from the folder containing ipwadhc.ado, ipwadh_covlists.ado, ipwadh_split.ado, ipwadh_mata.do,
* ipwadh_eval.ado (v0.2.0) and the folders data/S0_N2000.dta ... data/S7_N2000.dta (flat files, S3 not needed).   NOT YET RUN.
* For each scenario (N = 2000, seed 20261005; S3 was verified earlier): weights -> exact split -> truncation + diagnostics -> stcox, then one banded split (end),
* compared with the independent Python reference (ref_pipeline.py: Breslow, robust cluster SE x sqrt(G/(G-1)) = Stata vce(cluster id)).
* The checks do NOT stop at the first difference: each prints PASS/FAIL and the totals are listed at the end -- send me the whole log.
* RAM: exact split 6-12M rows (about 0.7 GB); every stcox takes about 1-2 minutes; expect 30-60 minutes in total (S0, S1, S5 are the largest).
* Tolerances: counts exact, coefficients/SEs 2e-4 absolute, weights/diagnostics 1e-4 relative (S4 and S6 have extreme weights: if only those fail, the
* Poisson optimizer is the first suspect -- compare the printed rate-model results).
* S8 (missing covariates, needs mi impute chained) is NOT included; S7 = whole-day times (ties).
clear all
set more off
global nfail = 0
global nchk = 0
local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead)"
local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
local sp "id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu)"

capture program drop chk
program define chk
    args label value expected tol mode
    * mode: abs (default) or rel
    if "`mode'" == "" local mode "abs"
    if "`mode'" == "abs" local dif = abs(`value' - `expected')
    else local dif = reldif(`value', `expected')
    global nchk = $nchk + 1
    if `dif' < `tol' {
        di as txt "  PASS  `label': " %14.8f `value' "  (python " %14.8f `expected' ")"
    }
    else {
        global nfail = $nfail + 1
        di as err "  FAIL  `label': " %14.8f `value' "  (python " %14.8f `expected' ", `mode' diff " %10.3e `dif' ")"
    }
end

capture program drop run_scen
program define run_scen
    args s nex cex sex lo hi ctr str n pt mean p50 mx ess mb nb cb sb
    local base "id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead)"
    local lists "stcoxvars(c.age i.male) adhpredictors(i.tsi i.abx c.biomarker)"
    local sp "id(id) time0(t0) time1(t1) adherence(adh) failure(inf) lastonly(dead ltfu)"
    di as res _n "=============== scenario `s' ==============="
    use "data/`s'_N2000.dta", clear
    ipwadhc, `base' `lists' parts replace
    ipwadh_split, `sp'
    chk "`s' exact split rows" _N `nex' 0.5
    ipwadh_eval, time0(t0) time1(t1) adherence(adh) arm(arm) truncate(1 99) replace
    scalar lo_s = r(trunc_lo)
    scalar hi_s = r(trunc_hi)
    local d_n = r(n_used)
    local d_pt = r(pt)
    local d_mean = r(mean)
    local d_p50 = r(p50)
    local d_max = r(max)
    local d_ess = r(ess)
    chk "`s' trunc lower bound" lo_s `lo' 1e-4 rel
    chk "`s' trunc upper bound" hi_s `hi' 1e-4 rel
    chk "`s' adherent pieces" `d_n' `n' 0.5
    chk "`s' person-time" `d_pt' `pt' 1e-6 rel
    chk "`s' mean weight (person-time)" `d_mean' `mean' 1e-4 rel
    chk "`s' median weight" `d_p50' `p50' 1e-4 rel
    chk "`s' max weight" `d_max' `mx' 1e-4 rel
    chk "`s' ESS fraction" `d_ess' `ess' 1e-4 rel
    stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
    stcox i.arm, breslow vce(cluster id)
    chk "`s' IPW exact coef" _b[1.arm] `cex' 2e-4
    chk "`s' IPW exact SE" _se[1.arm] `sex' 2e-4
    stset t1 [pw=_ipw_pswt] if adh == 1, time0(t0) failure(inf)
    stcox i.arm, breslow vce(cluster id)
    chk "`s' truncated coef" _b[1.arm] `ctr' 2e-4
    chk "`s' truncated SE" _se[1.arm] `str' 2e-4
    * banded split
    use "data/`s'_N2000.dta", clear
    quietly ipwadhc, `base' `lists' parts replace
    ipwadh_split, `sp' bandpoints(`mb')
    chk "`s' bandpoints(`mb') rows" _N `nb' 0.5
    ipwadh_eval, time0(t0) time1(t1) evalat(end) adherence(adh) replace
    stset t1 [pw=_ipw_psw] if adh == 1, time0(t0) failure(inf)
    stcox i.arm, breslow vce(cluster id)
    chk "`s' bandpoints(`mb') end coef" _b[1.arm] `cb' 2e-4
    chk "`s' bandpoints(`mb') end SE" _se[1.arm] `sb' 2e-4
end

timer clear
timer on 1
    run_scen S0 10964132 0.03476652 0.04283297 0.8496836434 1.1491800016 0.03502686 0.04280713 10943916 1333752.0916 0.99918306 0.99841817 1.37038149 0.99737255 1000 116198 0.03480643 0.04282800
    run_scen S1 9039567 -0.43851831 0.04588069 0.8471411934 1.1295499689 -0.43880671 0.04590851 9020154 1333083.0270 1.00029341 1.00188037 1.28778049 0.99771330 1000 109882 -0.43854840 0.04589390
    run_scen S2 7462105 -0.36414198 0.06242195 0.0131458164 5.5703857093 -0.36319849 0.05480660 7436459 1251093.2731 0.98512679 0.82013456 33.95976433 0.41633617 1000 114013 -0.36385726 0.06166019
    run_scen S4 6383628 -0.57684077 0.10592726 0.0014172656 10.2812499045 -0.41125656 0.06196133 6357032 1203783.4388 1.19172468 0.69108759 330.80474332 0.06013359 1000 110596 -0.58065909 0.10505238
    run_scen S5 11750964 -0.41965937 0.04847858 0.1876493569 1.9520892508 -0.42441663 0.04599226 11747501 1539560.6795 0.99779179 0.98919404 5.64174449 0.90968385 1000 61318 -0.41987168 0.04835959
    run_scen S6 6373930 -0.54095607 0.14180850 0.0000000251 2.9229649113 -0.40126334 0.05690487 6342566 1185265.7657 0.54442933 0.32650860 220.21223866 0.05728731 1000 115180 -0.54765537 0.13706553
    run_scen S7 1246343 -0.39254242 0.05335751 0.0163557876 5.3671703149 -0.39084396 0.04989262 1220729 1235258.0000 0.97907889 0.81914049 36.49706610 0.40645949 100 80895 -0.39442612 0.05236163

timer off 1
timer list 1
di as res _n "test_ipwadh_scenarios_N2000.do: $nchk checks, $nfail failed"
