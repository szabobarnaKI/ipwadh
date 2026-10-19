*! ipwadh_mi 0.2.1  IP weights for non-adherence on MULTIPLY IMPUTED data; the exact-split data are left in memory (a result frame) as an mi dataset (default FLONG) that mi estimate: stcox AND mi estimate: mean can use
*! NOT YET RUN IN STATA.  Design: PLAN_mi_splitdata.md; pieces tested in test_mi_stset_weights*.do and test_mi_build_smoke.do (wide / flong sources).
*!
*! WHAT IT DOES
*!   1. work on copies (frames): the user's mi data are never changed.  Every row gets a key _ipw_row (= _mi_id; = _n in wide data) so that the rows of every imputation can be put in the
*!      same order (mi extract m keeps the order for wide and flong source data, NOT for mlong).
*!   2. the m=0 data give the STRUCTURE (id, times, adherence, failure, death, ltfu, strata, arm: none of them may be imputed or passive).  ipwadh_split cuts the rows ONCE (exact split at the failure times of the adherent rows).
*!   3. for every imputation m = 1..M:  mi extract m -> assert the same structure as m=0 -> ipwadhc (weight models from THIS imputation's covariates, parts kept) -> the parts are mapped to the pieces by _ipw_orig
*!      -> ipwadh_eval (piece weights, truncation, diagnostics) -> psw, pswt and the person-time weights ptw = psw * (time1 - time0), ptwt = pswt * (time1 - time0) -> stored in the m-th imputation of the result.
*!   4. v0.2: DEFAULT result = FLONG mi data (mistyle(flong)): the weights are UNREGISTERED "super-varying" variables (m=0 value 1; mi update/mi stset/mi estimate accept them; measured in test_mi_stset_weights2.do C6a/C6b) so that
*!      mi stset [pw=_ipw_psw] ... and  mi estimate: stcox ...  work natively, and  mi estimate: mean x [pw=_ipw_ptw], over(arm)  too.  mistyle(wide) = the cheaper WIDE variant (weights registered imputed; tables only, no mi stset).
*!      (the text below describes the wide variant)  result frame: WIDE mi data (mi set wide).  Registered IMPUTED: the imputed and passive covariates of the original data (passive ones are copied, not re-derived; registered passive),
*!      _ipw_psw _ipw_ptw (and _ipw_pswt _ipw_ptwt with truncation); the m=0 value of the weights is missing (they cannot be computed from incomplete rows) so that  mi update  keeps the m>0 values.
*!      Everything else is registered regular.  No mi stset (mi stset refuses imputed weights; see PLAN section 2).
*!
*! USE
*!   weighted baseline / characteristics table, person-time weighted, adherent pieces:
*!       frame change <name>
*!       mi estimate: mean x1 x2 [pw=_ipw_ptw] if <adherence> == 1, over(<arm>)
*!       mi estimate: proportion g [pw=_ipw_ptw] if <adherence> == 1, over(<arm>)
*!   (_ipw_ptwt for the truncated weights).  Hazard ratios: loop over  mi extract m  -> stset t1 [pw=_ipw_psw] if adh==1, time0(t0) failure(..) -> stcox ...  -> Rubin's rules (see test_ipwadh_S8_mi.do).
*!
*! SYNTAX
*!   ipwadh_mi [if] [in], ID(varname) TIME0(varname) TIME1(varname) ADHerence(varname) ARM(varname) FAILure(varname) FRame(name [, replace])
*!        [ DEAD(varname) LTFU(varname) TERMinal(censor|none) STCOXvars(string) ADHPredictors(string) CENSORPredictors(string) BYARM POISsonopts(string)
*!          STRata(varlist) BANDpoints(# | auto) MAXGB(#) FORCE TRUNCate(# #) NOTRUNcate EVALat(end|mid|start) ESSCaution(#) ESSSevere(#) MAXWarn(#) CHange VERBose
*!          MIStyle(flong|wide) MIStset(psw|pswt) KEEP(varlist) ]
*!   The options are those of ipwadh (same meaning).  if / in are evaluated on the m=0 data (do not use imputed variables).  change = switch to the result frame at the end.
*!   Memory: the mi data are copied once per imputation (frame copy) in addition to the result.
*! RESULT characteristics (_dta[ipwadhmi_*] of the result frame): id time0 time1 adherence arm failure weight prefix M style of the source; no results are returned in r().

program define ipwadh_mi
    version 16
    local curframe = c(frame)
    capture noisily ipwadh_mi_run `0'
    local rc = _rc
    capture frame change `curframe'
    foreach f in _ipwmi_Fc _ipwmi_F0 _ipwmi_Fm {
        capture frame drop `f'
    }
    foreach mm in IPWMI_IDX IPWMI_SEL IPWMI_S0 IPWMI_SID IPWMI_PARTS IPWMI_VALS IPWMI_WT IPWMI_R {
        capture mata: mata drop `mm'
    }
    capture matrix drop IPWMI_DIAG
    if `rc' & "$IPWMI_FS" != "" {
        capture frame drop $IPWMI_FS
    }
    local fs "$IPWMI_FS"
    global IPWMI_FS ""
    if `rc' == 0 & "$IPWMI_CHANGE" == "1" {
        frame change `fs'
    }
    global IPWMI_CHANGE ""
    exit `rc'
end

program define ipwadh_mi_run
    version 16
    syntax [if] [in], ID(varname) TIME0(varname) TIME1(varname) ADHerence(varname) ARM(varname) FAILure(varname) FRame(string) ///
        [ DEAD(varname) LTFU(varname) TERMinal(string) STCOXvars(string) ADHPredictors(string) CENSORPredictors(string) BYARM POISsonopts(string) ///
          STRata(varlist) BANDpoints(string) MAXGB(real 8) FORCE ///
          TRUNCate(numlist min=2 max=2 >=0 <=100) NOTRUNcate EVALat(string) ESSCaution(real 20) ESSSevere(real 5) MAXWarn(real 100) ///
          CHange VERBose MIStyle(string) MIStset(string) KEEP(varlist) ]

    global IPWMI_FS ""
    global IPWMI_CHANGE ""
    local cur = c(frame)

    * ----- is it mi data? ---------------------------------------------------------------------------------------------------------------
    if "`_dta[_mi_marker]'" != "_mi_ds_1" {
        di as err "the data are not mi set. ipwadh_mi is for multiply imputed data (use ipwadh, keepdata for ordinary data)."
        exit 119
    }
    local style "`_dta[_mi_style]'"
    if !inlist("`style'", "wide", "mlong", "flong") {
        di as err "mi style `style' is not supported (use mi convert to wide, mlong or flong)"
        exit 198
    }
    local M = real("`_dta[_mi_M]'")
    if missing(`M') | `M' < 2 {
        di as err "at least 2 imputations are needed (M = `_dta[_mi_M]')"
        exit 198
    }
    local ivars "`_dta[_mi_ivars]'"
    local pvars "`_dta[_mi_pvars]'"
    local cov "`ivars' `pvars'"
    local cov = trim("`cov'")
    foreach v of local cov {
        capture confirm numeric variable `v'
        if _rc {
            di as err "imputed / passive variable `v' is a string variable: not supported"
            exit 198
        }
    }

    * ----- the structure may not vary across m ----------------------------------------------------------------------------------------------
    local struct "`id' `time0' `time1' `adherence' `failure' `dead' `ltfu' `arm' `strata'"
    foreach v of local struct {
        foreach w of local cov {
            if "`v'" == "`w'" {
                di as err "`v' is registered imputed or passive. The structure of the data (id, times, adherence, failure, death, ltfu, strata, arm) must not vary across imputations."
                di as err "Use the per-imputation loop with the full pipeline inside each imputation (see test_ipwadh_S8_mi.do)."
                exit 198
            }
        }
    }
    foreach w of local cov {
        if regexm(`"`if'"', "(^|[^A-Za-z0-9_])`w'([^A-Za-z0-9_]|$)") {
            di as err "the if condition uses the imputed / passive variable `w'. if / in are evaluated on the m=0 data: do not use imputed variables."
            exit 198
        }
    }

    * ----- the result frame ------------------------------------------------------------------------------------------------------------------
    gettoken fname frest : frame, parse(", ")
    local fname = trim("`fname'")
    local frest = trim(subinstr("`frest'", ",", " ", .))
    if "`frest'" != "" & "`frest'" != "replace" {
        di as err "frame(name [, replace]): `frest' not allowed"
        exit 198
    }
    confirm name `fname'
    capture confirm frame `fname'
    if !_rc {
        if "`frest'" != "replace" {
            di as err "frame `fname' already exists; use frame(`fname', replace)"
            exit 110
        }
        if "`fname'" == "`cur'" {
            di as err "the result frame cannot be the current frame"
            exit 198
        }
        frame drop `fname'
    }

    * ----- options (as in ipwadh) -----------------------------------------------------------------------------------------------------------
    if "`truncate'" != "" & "`notruncate'" != "" {
        di as err "truncate() and notruncate cannot be combined"
        exit 198
    }
    if "`truncate'" == "" & "`notruncate'" == "" local truncate "1 99"
    if "`evalat'" == "" local evalat end
    if "`mistyle'" == "" local mistyle flong
    if !inlist("`mistyle'", "wide", "flong") {
        di as err "mistyle() must be flong (default) or wide"
        exit 198
    }
    if "`mistset'" != "" {
        if !inlist("`mistset'", "psw", "pswt") {
            di as err "mistset() must be psw or pswt"
            exit 198
        }
        if "`mistyle'" == "wide" {
            di as err "mistset() needs mistyle(flong): mi stset refuses weights that differ by imputation in wide data"
            exit 198
        }
        if "`mistset'" == "pswt" & "`truncate'" == "" {
            di as err "mistset(pswt) needs truncated weights (not notruncate)"
            exit 198
        }
    }
    local nq = cond("`verbose'" == "", "quietly", "noisily")
    local dt = cond("`verbose'" == "", "", "detail")

    local opt1 "id(`id') time0(`time0') time1(`time1') adherence(`adherence') arm(`arm')"
    if "`dead'" != "" local opt1 "`opt1' dead(`dead')"
    if "`ltfu'" != "" local opt1 "`opt1' ltfu(`ltfu')"
    if "`terminal'" != "" local opt1 "`opt1' terminal(`terminal')"
    if "`stcoxvars'" != "" local opt1 `"`opt1' stcoxvars(`stcoxvars')"'
    if "`adhpredictors'" != "" local opt1 `"`opt1' adhpredictors(`adhpredictors')"'
    if "`censorpredictors'" != "" local opt1 `"`opt1' censorpredictors(`censorpredictors')"'
    if "`byarm'" != "" local opt1 "`opt1' byarm"
    if `"`poissonopts'"' != "" local opt1 `"`opt1' poissonopts(`poissonopts')"'

    local lastonly "`dead' `ltfu'"
    local opt2 "id(`id') time0(`time0') time1(`time1') adherence(`adherence') failure(`failure')"
    if "`strata'" != "" local opt2 "`opt2' strata(`strata')"
    if "`bandpoints'" != "" local opt2 "`opt2' bandpoints(`bandpoints')"
    if "`force'" != "" local opt2 "`opt2' force"
    local opt2 "`opt2' maxgb(`maxgb') `dt'"
    if trim("`lastonly'") != "" local opt2 "`opt2' lastonly(`lastonly')"

    local opt3 "time0(`time0') time1(`time1') adherence(`adherence') arm(`arm') evalat(`evalat') esscaution(`esscaution') esssevere(`esssevere') maxwarn(`maxwarn')"
    if "`truncate'" != "" local opt3 "`opt3' truncate(`truncate')"

    foreach s in orig piece row sel {
        capture confirm variable _ipw_`s'
        if !_rc {
            di as err "the data contain a variable _ipw_`s' (reserved by ipwadh)"
            exit 110
        }
    }

    * ----- 1. copy with a row key ------------------------------------------------------------------------------------------------------------
    di
    di as txt "ipwadh_mi: on-treatment (per-protocol) weights for non-adherence on mi data (" as res "`style'" as txt ", M = " as res `M' as txt ")"
    frame copy `cur' _ipwmi_Fc
    frame _ipwmi_Fc {
        if "`style'" == "wide" quietly gen long _ipw_row = _n
        else quietly gen long _ipw_row = _mi_id
    }

    * ----- 2. structure from m=0, ONE exact split ------------------------------------------------------------------------------------------
    frame copy _ipwmi_Fc _ipwmi_F0
    frame _ipwmi_F0 {
        quietly mi extract 0, clear
        sort _ipw_row
        quietly count
        local N0 = r(N)
        marksample touse, novarlist
        quietly gen byte _ipw_sel = `touse'
        quietly count if _ipw_sel == 1
        if r(N) == 0 error 2000
        local numstruct ""
        foreach v in `id' `time0' `time1' `adherence' `failure' `dead' `ltfu' `arm' `strata' {
            capture confirm numeric variable `v'
            if !_rc & strpos(" `numstruct' ", " `v' ") == 0 local numstruct `numstruct' `v'
        }
        local numstruct `numstruct' _ipw_row
        local strid = 0
        capture confirm string variable `id'
        if !_rc local strid = 1
        mata: IPWMI_S0 = st_data(., "`numstruct'")
        if `strid' mata: IPWMI_SID = st_sdata(., "`id'")
        mata: IPWMI_SEL = st_data(., "_ipw_sel")
        `nq' ipwadh_split if _ipw_sel == 1, `opt2'
        local n_pieces = _N
        mata: IPWMI_IDX = st_data(., "_ipw_orig")
    }

    * ----- 3. the result frame (wide mi data) ---------------------------------------------------------------------------------------------
    frame copy _ipwmi_F0 `fname'
    global IPWMI_FS "`fname'"
    frame `fname' {
        quietly gen double _ipw_psw = .
        quietly gen double _ipw_ptw = .
        label variable _ipw_psw "IPW at piece end, untruncated (imputed: differs by m)"
        label variable _ipw_ptw "person-time weight: _ipw_psw * (time1 - time0) (imputed)"
        local wvars "_ipw_psw _ipw_ptw"
        if "`truncate'" != "" {
            quietly gen double _ipw_pswt = .
            quietly gen double _ipw_ptwt = .
            label variable _ipw_pswt "IPW at piece end, truncated p`: word 1 of `truncate''/p`: word 2 of `truncate'' (imputed)"
            label variable _ipw_ptwt "person-time weight from the truncated weight (imputed)"
            local wvars "_ipw_psw _ipw_pswt _ipw_ptw _ipw_ptwt"
        }
        if "`keep'" != "" {
            local retain "`id' `time0' `time1' `adherence' `failure' `dead' `ltfu' `arm' `strata' `keep' `cov' `wvars' _ipw_row _ipw_orig _ipw_piece _ipw_sel"
            ds
            foreach v in `r(varlist)' {
                if strpos(" `retain' ", " `v' ") == 0 quietly drop `v'
            }
        }
        ds
        local allv "`r(varlist)'"
        local skipw = cond("`mistyle'" == "wide", "`wvars'", "")
        local regs ""
        foreach v of local allv {
            if strpos(" `cov' `skipw' ", " `v' ") == 0 local regs `regs' `v'
        }
        mi set wide
        if "`mistyle'" == "wide" mi register imputed `ivars' `wvars'
        else mi register imputed `ivars'
        if "`pvars'" != "" mi register passive `pvars'
        mi register regular `regs'
        mi set M = `M'
        if "`mistyle'" == "flong" {
            quietly mi convert flong, clear
            quietly mi unregister `wvars'
            quietly replace _ipw_psw = 1 if _mi_m == 0
            quietly replace _ipw_ptw = `time1' - `time0' if _mi_m == 0
            if "`truncate'" != "" {
                quietly replace _ipw_pswt = 1 if _mi_m == 0
                quietly replace _ipw_ptwt = `time1' - `time0' if _mi_m == 0
            }
            mata: st_local("okf", strofreal(rows(st_data(., "_mi_m")) == `n_pieces' * (`M' + 1) & all(st_data(., "_mi_m") :== vec(J(`n_pieces', 1, 1) * (0..`M'))) & all(st_data(., "_mi_id") :== vec((1::`n_pieces') * J(1, `M' + 1, 1)))))
            if `okf' != 1 {
                di as err "internal: the flong result has an unexpected row layout (_mi_m / _mi_id); please send this output"
                exit 459
            }
        }
    }

    * parts placeholders in the work frame (the pieces of the exact split)
    frame _ipwmi_F0 {
        foreach p in E0 J0 dA Ec0 dC {
            quietly gen double _ipw_`p' = .
        }
    }

    * ----- 4. every imputation ---------------------------------------------------------------------------------------------------------------
    matrix IPWMI_DIAG = J(`M', 4, .)
    local n_rows = .
    local n_changes = .
    local censoring ""
    forvalues m = 1/`M' {
        di as txt "  imputation " as res `m' as txt " of `M'" _c
        frame copy _ipwmi_Fc _ipwmi_Fm
        frame _ipwmi_Fm {
            quietly mi extract `m', clear
            sort _ipw_row
            quietly count
            local same = (r(N) == `N0')
            if `same' {
                mata: st_local("same", strofreal(all(st_data(., "`numstruct'") :== IPWMI_S0)))
            }
            if `same' & `strid' {
                mata: st_local("same", strofreal(all(st_sdata(., "`id'") :== IPWMI_SID)))
            }
            if `same' != 1 {
                di as err _n "the structure of imputation `m' differs from m=0 (rows, id, times, adherence, failure, death, ltfu, strata or arm): not allowed."
                exit 459
            }
            quietly gen byte _ipw_sel = .
            mata: st_store(., "_ipw_sel", IPWMI_SEL)
            `nq' ipwadhc if _ipw_sel == 1, `opt1' parts replace `dt'
            if `m' == 1 {
                local n_rows = r(N)
                local n_changes = r(n_changes)
                local censoring "`r(censoring)'"
            }
            sort _ipw_row
            mata: IPWMI_PARTS = st_data(., "_ipw_E0 _ipw_J0 _ipw_dA _ipw_Ec0 _ipw_dC")
            if "`cov'" != "" mata: IPWMI_VALS = st_data(., "`cov'")
        }
        frame drop _ipwmi_Fm

        frame _ipwmi_F0 {
            mata: st_store(., ("_ipw_E0", "_ipw_J0", "_ipw_dA", "_ipw_Ec0", "_ipw_dC"), IPWMI_PARTS[IPWMI_IDX, .])
            foreach dv in _ipw_psw _ipw_pswa _ipw_pswc _ipw_pswt _ipw_ptw _ipw_ptwt {
                capture drop `dv'
            }
            quietly ipwadh_eval if _ipw_sel == 1, `opt3' replace
            local d_ess = r(ess)
            local d_essmin = r(ess_min)
            local d_max = r(max)
            local d_warn = cond(missing(r(warn)), 0, r(warn))
            quietly gen double _ipw_ptw = _ipw_psw * (`time1' - `time0')
            local wlist "_ipw_psw _ipw_ptw"
            if "`truncate'" != "" {
                quietly gen double _ipw_ptwt = _ipw_pswt * (`time1' - `time0')
                local wlist "_ipw_psw _ipw_pswt _ipw_ptw _ipw_ptwt"
            }
            mata: IPWMI_WT = st_data(., "`wlist'")
        }
        matrix IPWMI_DIAG[`m', 1] = `d_ess'
        matrix IPWMI_DIAG[`m', 2] = `d_essmin'
        matrix IPWMI_DIAG[`m', 3] = `d_max'
        matrix IPWMI_DIAG[`m', 4] = `d_warn'

        frame `fname' {
            if "`mistyle'" == "wide" {
                local names ""
                foreach v of local wlist {
                    local names `names' _`m'_`v'
                }
                mata: st_store(., tokens("`names'"), IPWMI_WT)
                if "`cov'" != "" {
                    local names ""
                    foreach v of local cov {
                        local names `names' _`m'_`v'
                    }
                    mata: st_store(., tokens("`names'"), IPWMI_VALS[IPWMI_IDX, .])
                }
            }
            else {
                mata: IPWMI_R = (`m' * `n_pieces' + 1)::((`m' + 1) * `n_pieces')
                mata: st_store(IPWMI_R, tokens("`wlist'"), IPWMI_WT)
                if "`cov'" != "" mata: st_store(IPWMI_R, tokens("`cov'"), IPWMI_VALS[IPWMI_IDX, .])
            }
        }
        capture mata: mata drop IPWMI_PARTS
        capture mata: mata drop IPWMI_VALS
        capture mata: mata drop IPWMI_WT
        di as txt "   ESS " as res %4.1f 100*`d_ess' as txt "%, max weight " as res %7.2f `d_max'
    }

    * ----- 5. finish the result --------------------------------------------------------------------------------------------------------------
    frame `fname' {
        capture noisily mi update
        char define _dta[ipwadhmi] "1"
        char define _dta[ipwadhmi_id] "`id'"
        char define _dta[ipwadhmi_time0] "`time0'"
        char define _dta[ipwadhmi_time1] "`time1'"
        char define _dta[ipwadhmi_adherence] "`adherence'"
        char define _dta[ipwadhmi_arm] "`arm'"
        char define _dta[ipwadhmi_failure] "`failure'"
        char define _dta[ipwadhmi_M] "`M'"
        char define _dta[ipwadhmi_mistyle] "`mistyle'"
        if "`mistset'" != "" {
            local wv = cond("`mistset'" == "psw", "_ipw_psw", "_ipw_pswt")
            quietly mi stset `time1' [pw=`wv'] if `adherence' == 1, time0(`time0') failure(`failure')
        }
        char define _dta[ipwadhmi_sourcestyle] "`style'"
        char define _dta[ipwadhmi_truncate] "`truncate'"
        char define _dta[ipwadhmi_evalat] "`evalat'"
    }

    di
    di as txt "{hline 78}"
    di as txt "  imputation" _col(16) "ESS %" _col(28) "smallest arm ESS %" _col(50) "max weight" _col(64) "warning"
    forvalues m = 1/`M' {
        local w = IPWMI_DIAG[`m', 4]
        local wt = cond(`w' == 2, "SEVERE", cond(`w' == 1, "caution", ""))
        di as txt "  " %6.0f `m' _col(14) as res %7.1f 100*IPWMI_DIAG[`m', 1] _col(30) %7.1f 100*IPWMI_DIAG[`m', 2] _col(50) %9.2f IPWMI_DIAG[`m', 3] _col(64) as err "`wt'"
    }
    di as txt "{hline 78}"
    di as txt "  " as res %-9.0fc `n_rows' as txt " episodes, " as res %-8.0fc `n_changes' as txt " adherence changes, " as res %-9.0fc `n_pieces' as txt " pieces after the exact split (once, for all imputations)"
    di as txt "  result frame " as res "`fname'" as txt ": `mistyle' mi data, M = `M'" as txt "; switch:   frame change `fname'"
    if "`mistyle'" == "flong" {
        di as txt "  weights " as res "`wvars'" as txt " are UNREGISTERED super-varying variables (m=0 value: weight 1); " as res %-9.0fc `n_pieces' * (`M' + 1) as txt " rows stored"
        if "`mistset'" != "" di as txt "  mi stset done with [pw=_ipw_`mistset'] if `adherence' == 1"
        else di as txt "  mi stset `time1' [pw=_ipw_psw] if `adherence' == 1, time0(`time0') failure(`failure')      (or pswt)"
        di as txt "  mi estimate: stcox i.`arm' <covariates>, breslow vce(cluster `id')"
    }
    else di as txt "  weights " as res "`wvars'" as txt " are registered imputed (m=0 value missing); no mi stset / mi estimate: stcox in wide data (use mistyle(flong))"
    di as txt "  baseline table, person-time weighted, adherent pieces:"
    di as txt "    mi estimate: mean <x> [pw=_ipw_ptw] if `adherence' == 1, over(`arm')"
    di as txt "    mi estimate: proportion <g> [pw=_ipw_ptw] if `adherence' == 1, over(`arm')"
    if "`truncate'" != "" di as txt "    (truncated weights: _ipw_pswt / _ipw_ptwt)"
    di as txt "  note: if / in refer to the m=0 data; passive variables are copied, not re-derived. Memory: flong holds M+1 copies of every kept variable (keep() reduces)."

    if "`change'" != "" global IPWMI_CHANGE "1"
    global IPWMI_FS "`fname'"
end
