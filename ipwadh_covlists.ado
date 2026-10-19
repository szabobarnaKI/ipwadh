*! ipwadh_covlists 0.1.1  Plain-language covariate lists -> weight-model covariate lists (nesting guaranteed by construction)
*! 0.1.1: returned lists are trimmed (r(cdenvars) had a trailing blank when censorpredictors() was empty).
*! NOT YET FULLY RUN IN STATA (written in a container without Stata). Run test_ipwadh_covlists.do first.
*!
*! The user supplies TWO lists, no statistical jargon needed:
*!   stcoxvars(...)      the covariates you will put in your stcox / LWYY model, besides the treatment arm
*!                       (factor-variable notation allowed, e.g. c.age i.male i.arm#c.age)
*!   adhpredictors(...)  further variables that influence whether patients stay on / stop / restart treatment
*!                       and may also influence infections (time since infection, antibiotics, biomarker, ...)
*!   censorpredictors(...)  OPTIONAL: further predictors of death / loss to follow-up only
*!
*! WHAT IT DOES (so that the weights cannot be biased by a wrong list structure)
*!   - model for "reference" weights   = arm + current adherence + stcoxvars
*!   - model for "full" weights        = the reference model + adhpredictors   (always nested; stcoxvars are never dropped)
*!   - a variable already in stcoxvars is not repeated in adhpredictors (matched on the expanded variable names, so
*!     age, c.age and c.age#c.1 style duplicates are caught); partial overlaps (e.g. an interaction) are kept and flagged
*!   - the arm main effect is added automatically (any i.arm / c.arm / ib1.arm typed by the user is ignored); with byarm the
*!     models are fitted separately per arm, so arm and arm-interaction terms are ignored
*!   - missing values in any listed variable are reported by name (they would silently drop rows from the weight models)
*!   Returns for ipwadhc: r(numvars) r(denvars) r(cnumvars) r(cdenvars), and r(stcox_hint) = what to type in stcox.
*!   Bias note: this guarantees the structure, not the content: a confounder missing from adhpredictors is not repaired.

program define ipwadh_covlists, rclass
    version 15
    syntax [if] [in], ARM(varname) ADHerence(varname) [ STCOXvars(string) ADHPredictors(string) ///
        CENSORPredictors(string) BYARM ALLOWMissing QUIet ]

    marksample touse, novarlist
    local byflag = ("`byarm'" != "")

    local seen ""
    local numvars ""
    local denvars ""
    local cextra ""
    local n_arm ""
    local n_dup ""
    local n_part ""

    * ---- stcox covariates: go into the reference AND the full model ----
    foreach tok of local stcoxvars {
        _ipwadh_cl_term "`tok'" "`arm'" `byflag' "`seen'" "`touse'"
        if "`r(action)'" == "add" {
            local numvars `numvars' `tok'
            local seen `seen' `r(expanded)'
            if "`r(partial)'" == "1" local n_part `n_part' `tok'
        }
        else if "`r(action)'" == "skiparm" local n_arm `n_arm' `tok'
        else local n_dup `n_dup' `tok'
    }
    * ---- adherence predictors: full model only ----
    foreach tok of local adhpredictors {
        _ipwadh_cl_term "`tok'" "`arm'" `byflag' "`seen'" "`touse'"
        if "`r(action)'" == "add" {
            local denvars `denvars' `tok'
            local seen `seen' `r(expanded)'
            if "`r(partial)'" == "1" local n_part `n_part' `tok'
        }
        else if "`r(action)'" == "skiparm" local n_arm `n_arm' `tok'
        else local n_dup `n_dup' `tok'
    }
    * ---- extra death / loss predictors: full censoring model only ----
    foreach tok of local censorpredictors {
        _ipwadh_cl_term "`tok'" "`arm'" `byflag' "`seen'" "`touse'"
        if "`r(action)'" == "add" {
            local cextra `cextra' `tok'
            local seen `seen' `r(expanded)'
        }
        else if "`r(action)'" == "skiparm" local n_arm `n_arm' `tok'
        else local n_dup `n_dup' `tok'
    }

    * ---- missing values in any variable that will be used ----
    local alltok `numvars' `denvars' `cextra'
    local bad ""
    if "`alltok'" != "" {
        quietly fvrevar `alltok' if `touse', list
        local basevars `r(varlist)'
        local basevars : list uniq basevars
        foreach v of local basevars {
            quietly count if `touse' & missing(`v')
            if r(N) > 0 local bad `bad' `v'(`r(N)')
        }
    }
    if "`bad'" != "" & "`allowmissing'" == "" {
        di as err "missing values (number of rows in parentheses) in: `bad'"
        di as err "Rows with missing covariates would silently drop out of the weight models. Impute or recode them first,"
        di as err "or specify allowmissing to accept the loss of those rows (not recommended)."
        exit 416
    }

    if "`quiet'" == "" {
        local armtxt "treatment arm"
        if `byflag' local armtxt "(separate models for each arm)"
        di
        di as txt "Covariates used for the treatment-adherence and censoring weights"
        di as txt "{hline 72}"
        di as txt "  always, automatically:  `armtxt', current adherence state"
        di as txt "  your stcox covariates:  " as res cond("`numvars'" == "", "(none)", "`numvars'")
        di as txt "  added to explain why patients stop/restart treatment (kept out of the comparison):"
        di as txt "                          " as res cond("`denvars'" == "", "(none)", "`denvars'")
        if "`cextra'" != "" di as txt "  extra death/loss predictors: " as res "`cextra'"
        di as txt "{hline 72}"
        di as txt "  In your stcox, include the arm and exactly your stcox covariates:  " as res "i.arm `numvars'"
        if "`n_arm'" != "" di as txt "  ignored (arm is handled automatically): `n_arm'"
        if "`n_dup'" != "" di as txt "  not repeated (already in the model): `n_dup'"
        if "`n_part'" != "" di as txt "  note: partly overlaps with an earlier term, kept: `n_part'"
        if "`bad'" != "" di as err "  WARNING allowmissing: rows lost because of missing values in: `bad'"
    }

    * unquoted -local- assignments strip leading/trailing blanks
    local cnumvars `numvars'
    local cdenvars `denvars' `cextra'
    local hint i.arm `numvars'
    return local numvars "`numvars'"
    return local denvars "`denvars'"
    return local cnumvars "`cnumvars'"
    return local cdenvars "`cdenvars'"
    return local stcox_hint "`hint'"
    return local ignored_arm "`n_arm'"
    return local ignored_dup "`n_dup'"
end

* classify one term: add / skiparm / skipdup (all expanded names already present)
program define _ipwadh_cl_term, rclass
    args tok arm byflag seen touse
    if regexm("`tok'", "^((i|c|ib[0-9]*|ibn)\.)?`arm'$") {
        return local action "skiparm"
        exit
    }
    if `byflag' & regexm("`tok'", "(^|[.#])`arm'($|#)") {
        return local action "skiparm"
        exit
    }
    capture fvexpand `tok' if `touse'
    if _rc {
        di as err "could not interpret '`tok'' as a variable or factor-variable term (check spelling and that the variable exists)"
        exit 198
    }
    local ex `r(varlist)'
    local new : list ex - seen
    if "`new'" == "" {
        return local action "skipdup"
        exit
    }
    return local action "add"
    return local expanded "`ex'"
    local old : list ex & seen
    return local partial = ("`old'" != "")
end
