* test_ipwadhc_variants.do -- run in Stata 19 from the folder containing ipwadhc.ado (v0.2.2), ipwadh_covlists.ado, expected_ipwadhc_variants_S3_N200.csv
* and data/S3_N200/episodes.dta, data/S3_N2000/episodes.dta.   NOT YET RUN.
* Compares ipwadhc options  byarm  and  truncate(1 99)  (and the default models) with the independent Python reference (ref_ipwadhc.py with byarm=True,
* truncation by the Stata _pctile rule on ROW-weighted combined weights, ref_ipwadhc_variants.py):
*   S3_N200  : row by row (6,754 usable rows) swa swc sw swt, default and byarm, tolerance 1e-5 (relative); plus summary statistics
*   S3_N2000 : summary statistics (mean, sd, min, max of swa swc sw swt), truncation bounds, row and change counts, tolerance 1e-4 (relative)
* N200 byarm with Stata's DEFAULT Poisson convergence differed from Python by up to 5e-5 in the death weights (arm 1 has 11 deaths, none among antibiotic rows:
* nearly flat likelihood); the loose tolerance 2e-4 is used for that case, and the last section re-runs both modes with TIGHT Poisson tolerances
* (poissonopts(), ipwadhc 0.2.2) which must agree with Python to 1e-6. Watch for 'note:'/'warning:' lines printed by ipwadhc (perfect prediction / no convergence).
* Checks do NOT stop at the first difference: PASS/FAIL per check, totals at the end. Send me the whole log.
* Also re-run test_ipwadhc.do once (ipwadhc 0.2.1 only changed the closing advice text and added r(trunc_lo)/r(trunc_hi)); it should still pass.
clear all
set more off
global nfail = 0
global nchk = 0

capture program drop chk
program define chk
    args label value expected tol mode
    * mode: abs (default), rel (Stata reldif), rel2 (|value - expected| / |expected|)
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

capture program drop cmp_all
program define cmp_all
    args tag mode lo hi n nch popts
    local tol = cond("`tag'" == "N200", 1e-5, 1e-4)
    if "`tag'" == "N200" & "`mode'" == "byarm" local tol = 2e-4
    local po ""
    if `"`popts'"' != "" {
        local tol = 1e-6
        local po `"poissonopts(`popts')"'
    }
    di as res _n "=============== S3_`tag', `mode' `po' ==============="
    use "data/S3_`tag'/episodes.dta", clear
    local opt ""
    if "`mode'" == "byarm" local opt "byarm"
    ipwadhc, id(id) time0(t0) time1(t1) adherence(adh) arm(arm) dead(dead) stcoxvars(c.age i.male) ///
        adhpredictors(i.tsi i.abx c.biomarker) `opt' `po' truncate(1 99) replace
    local rN = r(N)
    local rch = r(n_changes)
    local wvar "`r(weightvar)'"
    scalar lo_s = r(trunc_lo)
    scalar hi_s = r(trunc_hi)
    if "`wvar'" != "_ipw_swt" {
        global nfail = $nfail + 1
        di as err "  FAIL  r(weightvar) is `wvar', expected _ipw_swt"
    }
    chk "`tag' `mode' usable rows" `rN' `n' 0.5
    chk "`tag' `mode' adherence changes" `rch' `nch' 0.5
    chk "`tag' `mode' truncation lower bound (p1)" lo_s `lo' `tol' rel2
    chk "`tag' `mode' truncation upper bound (p99)" hi_s `hi' `tol' rel2
    if "`tag'" == "N200" {
        preserve
        local pre "x_"
        if "`mode'" == "byarm" local pre "x_b_"
        merge 1:1 id rowno using "expected_variants_tmp.dta", keep(master match) generate(mg)
        quietly count if mg == 3
        chk "`tag' `mode' rows matched with the Python reference" r(N) 6754 0.5
        quietly count if (mg == 1) != (t1 == t0)
        chk "`tag' `mode' unmatched rows are exactly the zero-length rows" r(N) 0 0.5
        foreach v in swa swc sw swt {
            quietly gen double rd_`v' = abs(_ipw_`v' - `pre'`v') / abs(`pre'`v') if mg == 3
            quietly summarize rd_`v'
            local mx = r(max)
            chk "`tag' `mode' row-level max relative difference, `v'" `mx' 0 `tol'
        }
        restore
    }
end

capture program drop cmp_var
program define cmp_var
    args tag mode v m s mn mx
    local tol = cond("`tag'" == "N200", 1e-5, 1e-4)
    if "`tag'" == "N200" & "`mode'" == "byarm" local tol = 2e-4
    quietly summarize _ipw_`v'
    local a1 = r(mean)
    local a2 = r(sd)
    local a3 = r(min)
    local a4 = r(max)
    chk "`tag' `mode' `v' mean" `a1' `m'  `tol' rel2
    chk "`tag' `mode' `v' sd"   `a2' `s'  `tol' rel2
    chk "`tag' `mode' `v' min"  `a3' `mn' `tol' rel2
    chk "`tag' `mode' `v' max"  `a4' `mx' `tol' rel2
end

* expected row-level weights from the Python reference
import delimited using "expected_ipwadhc_variants_S3_N200.csv", clear
foreach v in swa swc sw swt b_swa b_swc b_sw b_swt {
    rename `v' x_`v'
}
assert _N == 6754
save "expected_variants_tmp.dta", replace

timer clear
timer on 1
    cmp_all N200 default 0.0058957589 7.0197614210 6754 2865
    cmp_var N200 default swa 1.0255625794 1.2873346997 0.0008512154 16.0634335948
    cmp_var N200 default swc 1.0004076291 0.0077668771 0.9537303486 1.0608626156
    cmp_var N200 default sw 1.0265847283 1.2920525105 0.0008517564 15.8237104706
    cmp_var N200 default swt 1.0059325200 1.1609787954 0.0058957589 7.0197614210
    cmp_all N200 byarm 0.0064953368 6.5164981556 6754 2865
    cmp_var N200 byarm swa 0.9867024345 1.2247033292 0.0004317998 19.7924390218
    cmp_var N200 byarm swc 1.0005485156 0.0108205634 0.9421154728 1.0614559915
    cmp_var N200 byarm sw 0.9872929025 1.2235442034 0.0004252665 19.0895871690
    cmp_var N200 byarm swt 0.9641970341 1.0684684158 0.0064953368 6.5164981556
    cmp_all N2000 default 0.0084212686 5.9910086298 69572 28775
    cmp_var N2000 default swa 0.9858883573 1.3416771578 0.0000510075 40.2516590300
    cmp_var N2000 default swc 1.0013858181 0.0173281324 0.8322327107 1.1425897825
    cmp_var N2000 default sw 0.9865747864 1.3551249286 0.0000534610 41.0805376600
    cmp_var N2000 default swt 0.9474051538 0.9953503630 0.0084212686 5.9910086298
    cmp_all N2000 byarm 0.0086324639 6.0099157027 69572 28775
    cmp_var N2000 byarm swa 0.9872029941 1.3616334197 0.0000472609 41.2693933287
    cmp_var N2000 byarm swc 1.0014027853 0.0182240170 0.8379806994 1.1786795560
    cmp_var N2000 byarm sw 0.9885875321 1.3824944713 0.0000487557 42.2298197094
    cmp_var N2000 byarm swt 0.9471794423 0.9965435053 0.0086324639 6.0099157027

* tight Poisson tolerances: the real test of the byarm death weights at N=200 (expected agreement with Python to about 1e-7)
    cmp_all N200 default 0.0058957589 7.0197614210 6754 2865 "tolerance(1e-12) ltolerance(1e-14) nrtolerance(1e-12)"
    cmp_all N200 byarm 0.0064953368 6.5164981556 6754 2865 "tolerance(1e-12) ltolerance(1e-14) nrtolerance(1e-12)"

erase "expected_variants_tmp.dta"
timer off 1
timer list 1
di as res _n "test_ipwadhc_variants.do: $nchk checks, $nfail failed"
