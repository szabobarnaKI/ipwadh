*! ipwadh 0.2.1  One-step wrapper: IP weights for non-adherence (and death / loss to follow-up) -> exact adherent-only split -> weights at failure times
*!               -> stset -> weighted Cox / LWYY model.  weightsonly = recommended route (then your own stcox); without it the command is a SHORTCUT that also fits the standard Cox model. Every step is also available as its own command (NOT YET RUN IN STATA)
*! 0.2.1: frame(name [, replace]) and change ALSO FOR ORDINARY DATA: a copy of the data is run with weightsonly in the new frame (your data are not touched; no model is fitted; the result frame is stset). NOT YET RUN IN STATA.
*! 0.2.0: MULTIPLY IMPUTED (mi set) data: with frame(name [, replace]) the call is handed to ipwadh_mi (result = mi dataset with the split data and weights in that frame; options MIStyle() MIStset() MIKeep() change). The non-mi path below is unchanged. NOT YET RUN IN STATA.
*! 0.1.1 passed test_ipwadh_wrapper.do in Stata 19 (32/33; the one failure was a wrong expected cluster count in the test).
*! 0.1.3: weightsonly ends with the Cox model to fit (arm + stcoxvars, strata()).
*! 0.1.2: compact output: one header line, the compact ipwadh_eval report with the weight warning, the Cox table (no iteration log, no stset echo); verbose = detail output of every step.
*!
*! WHAT IT DOES (in this order; each step is the stand-alone command of the same name)
*!   1. ipwadhc        rate models (Poisson) and the stabilized weights, with the `parts` needed for the evaluation inside rows
*!   2. ipwadh_split   exact adherent-only split at the failure times
*!   3. ipwadh_eval    weight of every piece (end of the piece), diagnostics, truncation (default p1/p99, person-time weighted), weight warning
*!   4. stset          t1 [pw=<weight>] if adherence==1, time0(t0) failure(<failure>)       (no id(), no exit())
*!   5. stcox          i.arm <stcoxvars>, breslow vce(cluster id)             DEFAULT: cluster-robust SE (weights treated as known)
*!                     with adjustse: stcox_ipwadh = SE that also accounts for the estimation of the weights (supplementary; differs by a few percent)
*!
*! The Cox model is always  i.<arm> + <stcoxvars>  (the rule that makes the weights valid for that model).  With strata(), the strata variables are used in the split
*! and in stcox (as strata(); they are dropped from the covariates of stcox but are kept in the weight models).
*!
*! DATA: by default the original (unsplit) data are RESTORED at the end; only the estimation results remain (e()).  Option keepdata keeps the split data, the
*! weights (_ipw_psw, _ipw_pswt) and the stset in memory (needed for predict, estat, other models, stcox_ipwadh afterwards, ipwadh_eval diagnostics).
*!
*! MULTIPLY IMPUTED DATA (0.2.0): if the data in memory are mi set, ipwadh needs  frame(name [, replace])  and calls ipwadh_mi with the same options (same meaning).  Your mi data are never changed.
*!   The result is an mi dataset in frame `name' (default FLONG): the exact-split data with the weights _ipw_psw _ipw_pswt (and the person-time weights _ipw_ptw _ipw_ptwt) as UNREGISTERED super-varying
*!   variables, so that  mi estimate: stcox  (after  mi stset ... [pw=_ipw_psw] / MIStset())  and  mi estimate: mean x [pw=_ipw_ptw] if adh==1, over(arm)  work directly.  No model is fitted by ipwadh itself.
*!   MIStyle(flong|wide)  MIStset(psw|pswt)  MIKeep(varlist)  change  = the options style / stset / keep / change of ipwadh_mi.  adjustse is not available; keepdata and weightsonly are implied; nohr and level are ignored.
*!   (the option is MIKeep() here because keep() would be ambiguous with keepdata).  MIStyle(), MIStset() and MIKeep() are an ERROR on ordinary data.
*!
*! RESULT IN A FRAME, ORDINARY DATA (0.2.1): frame(name [, replace]) copies the data to a new frame and runs  ipwadh ..., weightsonly  there: the split data, the weights and the stset stay in that frame, your data in the current frame are
*!   unchanged.  No model is fitted (adjustse is an error); change switches to the result frame.  Memory: the data exist twice.  An existing frame needs replace; the current frame cannot be the result frame; a failed call leaves no frame behind.
*!
*! SYNTAX
*!   ipwadh [if] [in], ID(varname) TIME0(varname) TIME1(varname) ADHerence(varname) ARM(varname) FAILure(varname)
*!        [ DEAD(varname) LTFU(varname) TERMinal(censor|none) STCOXvars(string) ADHPredictors(string) CENSORPredictors(string) BYARM POISsonopts(string)
*!          STRata(varlist) BANDpoints(# | auto) MAXGB(#) FORCE
*!          TRUNCate(# #) NOTRUNcate EVALat(end|mid|start) ESSCaution(#) ESSSevere(#) MAXWarn(#)
*!          ADJustse NOHR Level(#) KEEPdata WEIGHTSonly VERBose
*!          FRame(name [, replace]) MIStyle(flong|wide) MIStset(psw|pswt) MIKeep(varlist) CHange ]       (FRame, CHange: any data; MIStyle, MIStset, MIKeep: mi data only)
*!   truncate(1 99) is the default; notruncate uses the untruncated weights.  weightsonly stops after stset (implies keepdata): fit your own model.
*!   verbose shows the output of the individual steps.  Composite event (infection or death): define the event variable yourself, use it as failure(), terminal(none).
*! RETURNS (e-class, in addition to the stcox results): e(ipwadh_ess) e(ipwadh_ess_min) e(ipwadh_wmax) e(ipwadh_warn) e(ipwadh_ess_trunc) e(ipwadh_n_pieces) e(ipwadh_n_rows)
*!   e(ipwadh_weight) (name of the weight variable)  e(ipwadh_se) (cluster-robust | adjusted)       (mi data: nothing is returned; the result frame carries _dta[ipwadhmi_*] characteristics)

program define ipwadh, eclass
    version 15
    syntax [if] [in], ID(varname) TIME0(varname) TIME1(varname) ADHerence(varname) ARM(varname) FAILure(varname) ///
        [ DEAD(varname) LTFU(varname) TERMinal(string) STCOXvars(string) ADHPredictors(string) CENSORPredictors(string) BYARM POISsonopts(string) ///
          STRata(varlist) BANDpoints(string) MAXGB(real 8) FORCE ///
          TRUNCate(numlist min=2 max=2 >=0 <=100) NOTRUNcate EVALat(string) ESSCaution(real 20) ESSSevere(real 5) MAXWarn(real 100) ///
          ADJustse NOHR Level(cilevel) KEEPdata WEIGHTSonly VERBose ///
          FRame(string) MIStyle(string) MIStset(string) MIKeep(varlist) CHange ]

    * ----- result in a frame (0.2.0 / 0.2.1): mi data -> ipwadh_mi; ordinary data -> the weightsonly pipeline run in a COPY of the data (everything below this block is the unchanged ordinary path) ---
    local ismi = ("`_dta[_mi_marker]'" == "_mi_ds_1")
    if !`ismi' & ("`mistyle'" != "" | "`mistset'" != "" | "`mikeep'" != "") {
        di as err "mistyle(), mistset() and mikeep() are for multiply imputed (mi set) data"
        exit 198
    }
    if "`change'" != "" & "`frame'" == "" {
        di as err "change needs frame()"
        exit 198
    }
    if `ismi' & "`frame'" == "" {
        di as err "the data are mi set: give the result frame, frame(name [, replace]). Your mi data are not changed; the split data and the weights are returned as an mi dataset in that frame (see help ipwadh_mi)."
        exit 198
    }
    if "`frame'" != "" {
        if c(stata_version) < 16 {
            di as err "frame() needs Stata 16 or later"
            exit 198
        }
        if "`adjustse'" != "" {
            di as err "adjustse cannot be combined with frame(): no model is fitted, the result frame holds the split data and the weights (fit stcox or stcox_ipwadh there)"
            exit 198
        }
        * the options common to both routes
        local ropt `"id(`id') time0(`time0') time1(`time1') adherence(`adherence') arm(`arm') failure(`failure')"'
        foreach o in dead ltfu terminal stcoxvars adhpredictors censorpredictors poissonopts strata bandpoints evalat {
            if `"``o''"' != "" local ropt `"`ropt' `o'(``o'')"'
        }
        local ropt `"`ropt' maxgb(`maxgb') esscaution(`esscaution') esssevere(`esssevere') maxwarn(`maxwarn')"'
        if "`truncate'" != "" local ropt `"`ropt' truncate(`truncate')"'
        foreach o in notruncate byarm force verbose {
            if "``o''" != "" local ropt `"`ropt' `o'"'
        }
        if `ismi' {
            local mopt `"`ropt' frame(`frame')"'
            foreach o in mistyle mistset {
                if `"``o''"' != "" local mopt `"`mopt' `o'(``o'')"'
            }
            if "`change'" != "" local mopt `"`mopt' change"'
            if "`mikeep'" != "" local mopt `"`mopt' keep(`mikeep')"'
            if "`keepdata'`weightsonly'" != "" & "`verbose'" != "" di as txt "note: mi data: the result frame always holds the split data and the weights (keepdata / weightsonly are implied)"
            ipwadh_mi `if' `in', `mopt'
            exit
        }
        * ordinary data: copy the data to the new frame and run  ipwadh ..., weightsonly  there; the user's data are never touched
        gettoken fname frest : frame, parse(", ")
        local fname = trim("`fname'")
        local frest = trim(subinstr("`frest'", ",", " ", .))
        if "`frest'" != "" & "`frest'" != "replace" {
            di as err "frame(name [, replace]): `frest' not allowed"
            exit 198
        }
        confirm name `fname'
        local cur = c(frame)
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
        frame copy `cur' `fname'
        capture noisily frame `fname': ipwadh `if' `in', `ropt' weightsonly
        local rc = _rc
        if `rc' {
            capture frame drop `fname'
            exit `rc'
        }
        di as txt "result frame " as res "`fname'" as txt ": the split data, the weights and the stset; your data in frame " as res "`cur'" as txt " are unchanged.  Switch:   frame change `fname'"
        if "`keepdata'" != "" & "`verbose'" != "" di as txt "note: with frame() keepdata / weightsonly are implied and no model is fitted"
        if "`change'" != "" frame change `fname'
        exit
    }

    if "`truncate'" != "" & "`notruncate'" != "" {
        di as err "truncate() and notruncate cannot be combined"
        exit 198
    }
    if "`truncate'" == "" & "`notruncate'" == "" local truncate "1 99"
    if "`weightsonly'" != "" & "`adjustse'" != "" {
        di as err "weightsonly and adjustse cannot be combined"
        exit 198
    }
    if "`evalat'" == "" local evalat end
    if "`evalat'" != "end" & "`bandpoints'" == "" {
        di as txt "note: evalat(`evalat') is only needed with a banded split (bandpoints() > 1); the exact split uses evalat(end)"
    }
    if "`weightsonly'" != "" local keepdata keepdata
    foreach s in orig piece {
        capture confirm variable _ipw_`s'
        if !_rc {
            di as err "variable _ipw_`s' already exists: the data look already split. ipwadh needs the UNSPLIT episode data (use the individual commands to continue from a split)."
            exit 110
        }
    }
    capture confirm variable _ipw_psw
    if !_rc & "`verbose'" == "" di as txt "note: existing _ipw_* variables are replaced"
    local nq = cond("`verbose'" == "", "quietly", "noisily")

    marksample touse, novarlist

    * ----- the Cox covariates: arm + stcoxvars (strata variables are removed from the covariates) ----------------------------------------
    local coxcov ""
    foreach w of local stcoxvars {
        local base = subinstr(subinstr("`w'", "i.", "", .), "c.", "", .)
        local isstr 0
        foreach s of local strata {
            if "`base'" == "`s'" local isstr 1
        }
        if !`isstr' local coxcov `coxcov' `w'
    }
    local coxvars i.`arm' `coxcov'

    if "`keepdata'" == "" preserve

    * ----- 1. weights -------------------------------------------------------------------------------------------------------------------
    local opt1 "id(`id') time0(`time0') time1(`time1') adherence(`adherence') arm(`arm')"
    if "`dead'" != "" local opt1 "`opt1' dead(`dead')"
    if "`ltfu'" != "" local opt1 "`opt1' ltfu(`ltfu')"
    if "`terminal'" != "" local opt1 "`opt1' terminal(`terminal')"
    if "`stcoxvars'" != "" local opt1 `"`opt1' stcoxvars(`stcoxvars')"'
    if "`adhpredictors'" != "" local opt1 `"`opt1' adhpredictors(`adhpredictors')"'
    if "`censorpredictors'" != "" local opt1 `"`opt1' censorpredictors(`censorpredictors')"'
    if "`byarm'" != "" local opt1 "`opt1' byarm"
    if `"`poissonopts'"' != "" local opt1 `"`opt1' poissonopts(`poissonopts')"'
    local dt = cond("`verbose'" == "", "", "detail")
    `nq' ipwadhc if `touse', `opt1' parts replace `dt'
    local n_rows = r(N)
    local n_changes = r(n_changes)
    local censoring "`r(censoring)'"
    local terminal_used "`r(terminal)'"

    * ----- 2. split ---------------------------------------------------------------------------------------------------------------------
    local lastonly "`dead' `ltfu'"
    local opt2 "id(`id') time0(`time0') time1(`time1') adherence(`adherence') failure(`failure')"
    if "`strata'" != "" local opt2 "`opt2' strata(`strata')"
    if "`bandpoints'" != "" local opt2 "`opt2' bandpoints(`bandpoints')"
    if "`force'" != "" local opt2 "`opt2' force"
    local opt2 "`opt2' maxgb(`maxgb') `dt'"
    if trim("`lastonly'") != "" local opt2 "`opt2' lastonly(`lastonly')"
    `nq' ipwadh_split if `touse', `opt2'
    local n_pieces = _N
    if "`verbose'" == "" {
        di
        di as txt "ipwadh: on-treatment (per-protocol) analysis, weights for non-adherence" _c
        if "`censoring'" != "" & "`censoring'" != "nothing" di as txt " + censoring by `censoring'" _c
        di as txt ""
        di as txt "  " as res %-9.0fc `n_rows' as txt " episodes, " as res %-8.0fc `n_changes' as txt " adherence changes, " as res %-9.0fc `n_pieces' as txt " pieces after the exact split"
    }

    * ----- 3. weights at the pieces, diagnostics ----------------------------------------------------------------------------------------
    local opt3 "time0(`time0') time1(`time1') adherence(`adherence') arm(`arm') evalat(`evalat') esscaution(`esscaution') esssevere(`esssevere') maxwarn(`maxwarn')"
    if "`truncate'" != "" local opt3 "`opt3' truncate(`truncate')"
    if "`verbose'" != "" local opt3 "`opt3' detail"
    ipwadh_eval if `touse', `opt3' replace
    local ess = r(ess)
    local ess_min = r(ess_min)
    local wmax = r(max)
    local warn = cond(missing(r(warn)), 0, r(warn))
    local warn_max = cond(missing(r(warn_max)), 0, r(warn_max))
    local wv "`r(weightvar)'"
    local ess_t = .
    local share_lo = .
    local share_hi = .
    if "`truncate'" != "" {
        local ess_t = r(ess_trunc)
        local share_lo = r(share_lo)
        local share_hi = r(share_hi)
    }

    * ----- 4. stset ---------------------------------------------------------------------------------------------------------------------
    quietly stset `time1' [pw=`wv'] if `adherence' == 1 & `touse', time0(`time0') failure(`failure')

    * ----- 5. the Cox model -----------------------------------------------------------------------------------------------------------------
    if "`weightsonly'" == "" {
        local copt ""
        if "`strata'" != "" local copt "strata(`strata')"
        if "`adjustse'" != "" {
            stcox_ipwadh `coxvars', `copt' `nohr' level(`level') evalat(`evalat')
            local sekind "adjusted"
        }
        else {
            di
            stcox `coxvars', breslow vce(cluster `id') `copt' `nohr' level(`level') nolog noshow
            local sekind "cluster-robust"
        }
        ereturn scalar ipwadh_ess = `ess'
        ereturn scalar ipwadh_ess_min = `ess_min'
        ereturn scalar ipwadh_wmax = `wmax'
        ereturn scalar ipwadh_warn = `warn'
        ereturn scalar ipwadh_ess_trunc = `ess_t'
        ereturn scalar ipwadh_n_pieces = `n_pieces'
        ereturn scalar ipwadh_n_rows = `n_rows'
        ereturn local ipwadh_weight "`wv'"
        ereturn local ipwadh_se "`sekind'"
        ereturn local ipwadh_cmd "ipwadh `coxvars'"
    }
    else {
        di
        di as txt "weights ready: data are stset with [pw=`wv'] for the adherent pieces (split data left in memory)."
        local scov "`coxvars'"
        local sopt "breslow vce(cluster `id')"
        if "`strata'" != "" local sopt "`sopt' strata(`strata')"
        di as txt "  Cox model to fit:  stcox " as res "`scov', `sopt'"
        di as txt "  (arm plus exactly the stcoxvars(); adhpredictors() / censorpredictors() stay out)"
    }
    if "`keepdata'" != "" & "`weightsonly'" == "" di as txt "(keepdata: the split data, the weights and the stset are left in memory)"
end
