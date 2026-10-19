* test_mi_stset_weights.do  --  does  mi stset ... [pw=<varying weight>]  work?  (NOT YET RUN; written without access to Stata)
* Question: can the weight in  mi stset  be a registered IMPUTED or PASSIVE variable (values differ across m), as the split ipwadh mi data would need?
* Self-contained: synthetic data (200 persons x 2 episode rows, row-level imputed covariate x, 5 imputations), no files, no internet.  Runs in about a minute.  Stata 16+.
* For every case the do-file reports (a) whether  mi stset  was refused (rc) and (b) if not, whether  mi estimate: stcox  really uses the m-specific weights:
*   the pooled log HR from  mi estimate  is compared with the mean over m of the log HRs from  mi extract m -> stset [pw=] -> stcox  (reference, plain Stata).
* Styles tested: wide and flong (mlong cannot hold weights that vary on complete rows).
* Cases per style
*   C1  imputed weight, m=0 value missing        (the design of PLAN_mi_splitdata.md)
*   C1b imputed weight, m=0 value NOT missing    (does  mi update  overwrite the m>0 values?)
*   C2  passive weight  (mi passive: generate)
*   C3  regular weight  (control: must work; same in every m)
*   C4  workaround b:  passive copy -> mi unregister -> mi stset -> mi register passive   (Statalist svyset thread)
*   C5  workaround c:  regular placeholder -> mi stset -> mi register imputed -> overwrite m>0 values
clear all
set more off
version 16
global ncase = 0

* ---------------------------------------------------------------------------------------------------------------------------
capture program drop mkdata
program define mkdata
    args style
    clear
    set seed 4711
    set obs 200
    gen long id = _n
    gen byte arm = (runiform() < .5)
    gen double T = rexponential(8) * exp(-0.4 * arm)
    gen double time = min(T, 10)
    gen byte d = (T < 10)
    expand 2
    bysort id: gen byte row = _n
    gen double t0 = cond(row == 1, 0, time/2)
    gen double t1 = cond(row == 1, time/2, time)
    gen byte fail = (row == 2 & d == 1)
    gen double x = rnormal(0, 1) + 0.3 * arm
    replace x = . if runiform() < 0.25
    drop T time d
    sort id row
    mi set `style'
    mi register imputed x
    mi register regular id arm row t0 t1 fail
    quietly mi impute regress x = arm t1 fail, add(5) rseed(20261011)
end

* ---------------------------------------------------------------------------------------------------------------------------
* store one line for the final recap
capture program drop recap
program define recap
    args text
    global ncase = $ncase + 1
    global r$ncase `"`text'"'
end

* ---------------------------------------------------------------------------------------------------------------------------
* mi stset with the weight wvar; returns rc in r(rc)
capture program drop try_stset
program define try_stset, rclass
    args wvar
    capture noisily mi stset t1 [pw=`wvar'], time0(t0) failure(fail)
    return scalar rc = _rc
end

* ---------------------------------------------------------------------------------------------------------------------------
* data are mi stset already.  mi estimate vs per-m reference (weights wref in the extracted data)
capture program drop est_check
program define est_check
    args label wref
    capture noisily mi estimate, post: stcox i.arm, breslow vce(cluster id)
    local rc = _rc
    if `rc' {
        recap "`label': mi stset OK, but mi estimate: stcox FAILED, rc = `rc'"
        exit
    }
    local bmi = _b[1.arm]
    local sum = 0
    local bmin = .
    local bmax = .
    forvalues m = 1/5 {
        preserve
        quietly mi extract `m', clear
        quietly stset t1 [pw=`wref'], time0(t0) failure(fail)
        quietly stcox i.arm, breslow vce(cluster id)
        local b = _b[1.arm]
        local sum = `sum' + `b'
        if missing(`bmin') | `b' < `bmin' local bmin = `b'
        if missing(`bmax') | `b' > `bmax' local bmax = `b'
        restore
    }
    local bref = `sum' / 5
    local dif = abs(`bmi' - `bref')
    local spread = `bmax' - `bmin'
    local verdict "MISMATCH (mi estimate does NOT use the m-specific weights)"
    if `dif' < 1e-8 & `spread' > 1e-6 local verdict "OK: mi estimate uses the m-specific weights"
    if `dif' < 1e-8 & `spread' <= 1e-6 local verdict "agrees, but weights hardly differ across m (uninformative)"
    local line = "`label': mi stset OK;  mi estimate log HR " + string(`bmi', "%10.6f") + "  reference mean " + string(`bref', "%10.6f") + "  diff " + string(`dif', "%9.2e") + "  spread over m " + string(`spread', "%9.2e") + "  -> `verdict'"
    recap `"`line'"'
end

* ---------------------------------------------------------------------------------------------------------------------------
capture program drop show_m
program define show_m
    args wvar
    foreach m in 0 1 2 {
        quietly mi xeq `m': summarize `wvar'
        di as txt "    `wvar' in m=`m': N = " r(N) "  mean = " %9.5f r(mean) "  sd = " %9.5f r(sd)
    }
end

* ===========================================================================================================================
foreach style in wide flong {

    di as res _n "################  style `style'  ################"

    * ----- C1: imputed weight, m=0 missing -------------------------------------------------------------------------------
    di as res _n "--- `style' C1: imputed weight, m=0 missing"
    mkdata `style'
    gen double w_imp = .
    mi register imputed w_imp
    mi xeq 1/5: replace w_imp = exp(0.2 * x)
    capture noisily mi update
    show_m w_imp
    try_stset w_imp
    local rc_set = r(rc)
    if `rc_set' {
        recap "`style' C1  imputed weight, m=0 missing: mi stset REFUSED, rc = `rc_set'"
    }
    else {
        est_check "`style' C1  imputed weight, m=0 missing" w_imp
    }

    * ----- C1b: imputed weight, m=0 NOT missing --------------------------------------------------------------------------
    di as res _n "--- `style' C1b: imputed weight, m=0 value 1 (not missing)"
    mkdata `style'
    gen double w_imp2 = 1
    mi register imputed w_imp2
    mi xeq 1/5: replace w_imp2 = exp(0.2 * x)
    capture noisily mi update
    di as txt "  after mi update (are the m>0 values kept, or reset to the m=0 value 1?):"
    show_m w_imp2
    try_stset w_imp2
    local rc_set = r(rc)
    if `rc_set' {
        recap "`style' C1b imputed weight, m=0 = 1: mi stset REFUSED, rc = `rc_set'"
    }
    else {
        est_check "`style' C1b imputed weight, m=0 = 1" w_imp2
    }

    * ----- C2: passive weight --------------------------------------------------------------------------------------------
    di as res _n "--- `style' C2: passive weight"
    mkdata `style'
    mi passive: generate double w_pas = exp(0.2 * x)
    show_m w_pas
    try_stset w_pas
    local rc_set = r(rc)
    if `rc_set' {
        recap "`style' C2  passive weight: mi stset REFUSED, rc = `rc_set'"
    }
    else {
        est_check "`style' C2  passive weight" w_pas
    }

    * ----- C3: regular weight (control) ----------------------------------------------------------------------------------
    di as res _n "--- `style' C3: regular weight (control)"
    mkdata `style'
    gen double w_reg = 1 + 0.5 * arm
    mi register regular w_reg
    try_stset w_reg
    local rc_set = r(rc)
    if `rc_set' {
        recap "`style' C3  regular weight (control): mi stset REFUSED, rc = `rc_set'  <-- the test itself is broken"
    }
    else {
        est_check "`style' C3  regular weight (control)" w_reg
    }

    * ----- C4: workaround b: passive copy, unregister, stset, re-register passive -----------------------------------------
    di as res _n "--- `style' C4: passive copy -> mi unregister -> mi stset -> mi register passive"
    mkdata `style'
    mi passive: generate double w_pas = exp(0.2 * x)
    mi passive: generate double w_cp = w_pas
    capture noisily mi unregister w_cp
    local rc_unreg = _rc
    quietly mi xeq 1: count if abs(w_cp - w_pas) > 1e-9 & !missing(w_pas)
    local d1 = r(N)
    di as txt "  after unregister: rows in m=1 where the copy differs from the passive original: `d1'"
    try_stset w_cp
    local rc_set = r(rc)
    capture noisily mi register passive w_cp
    local rc_reg = _rc
    quietly mi xeq 1: count if abs(w_cp - w_pas) > 1e-9 & !missing(w_pas)
    local d2 = r(N)
    di as txt "  after re-register: rows in m=1 where the copy differs: `d2'"
    if `rc_set' recap "`style' C4  workaround b: mi stset REFUSED even for the unregistered copy, rc = `rc_set' (unregister rc `rc_unreg')"
    else {
        recap "`style' C4  workaround b: unregister rc `rc_unreg', stset rc 0, re-register rc `rc_reg'; copy differs from original in m=1 on `d1' rows after unregister, `d2' after re-register (0 = values survived)"
        est_check "`style' C4  workaround b (mi estimate)" w_pas
    }

    * ----- C5: workaround c: regular placeholder, stset, register imputed, overwrite m>0 -----------------------------------
    di as res _n "--- `style' C5: regular placeholder -> mi stset -> mi register imputed -> overwrite m>0 values"
    mkdata `style'
    gen double w_ph = 1
    mi register regular w_ph
    try_stset w_ph
    local rc_set = r(rc)
    if `rc_set' recap "`style' C5  workaround c: mi stset REFUSED for the regular placeholder, rc = `rc_set'  <-- the test itself is broken"
    else {
        capture noisily mi register imputed w_ph
        local rc_reg = _rc
        quietly mi xeq 0: replace w_ph = .
        quietly mi xeq 1/5: replace w_ph = exp(0.2 * x)
        capture noisily mi update
        di as txt "  after the overwrite and mi update:"
        show_m w_ph
        recap "`style' C5  workaround c: stset rc 0, register imputed rc `rc_reg'"
        est_check "`style' C5  workaround c (mi estimate)" w_ph
    }
}

* ===========================================================================================================================
di as res _n "=================================== RECAP ==================================="
forvalues i = 1/$ncase {
    di as txt "${r`i'}"
}
di as res _n "test_mi_stset_weights.do finished: send me the WHOLE log (the per-case blocks and this recap)."
