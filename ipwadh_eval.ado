*! ipwadh_eval 0.3.1  Attach IP weights to the pieces made by ipwadh_split; truncation, diagnostics and weight warnings (person-time weighted)
*! 0.1.0 PASSED test_ipwadh_eval.do in Stata 19.  0.2.0 (truncate(), arm(), person-time-weighted diagnostics): test_ipwadh_eval_trunc.do.
*! 0.3.0 adds the weight-instability warning (esscaution() esssevere() maxwarn() nowarn; returns r(warn) r(warn_max) r(ess_min) r(ess_trunc)): run test_ipwadh_eval_warn.do (NOT YET RUN).
*!
*! 0.3.1: compact default output (3 lines + warnings); option detail restores the full diagnostics table (by arm too) and the analysis hints. Numbers unchanged.
*!
*! WORKFLOW (order matters)
*!   1. ipwadhc, ... parts replace            on the UNSPLIT episodes: creates _ipw_E0 _ipw_J0 _ipw_dA _ipw_Ec0 _ipw_dC (+ row-end weights)
*!   2. ipwadh_split, ...                     cuts the rows; expand copies the parts onto every piece; creates _ipw_orig, _ipw_piece
*!   3. ipwadh_eval, time0(t0) time1(t1) ...  weight of every piece:
*!          log W(t) = E0 - dA (t - s) + J0 + Ec0 - dC (t - s)       s = start of the ORIGINAL row (time0 of piece 1), t = evalat point of the piece
*!   4. stset t1 [pw=_ipw_psw] if adh==1, time0(t0) failure(inf)     (NO id() and no exit(): Stata refuses time-varying weights with id(), r(459);
*!      every piece is its own record, the failure indicator is on the last piece only; one failure per record, so exit() is not needed)
*!      stcox i.arm, breslow vce(cluster id)
*!
*! evalat(end)    t = piece end (default). Exact for the exact split (bandpoints(1)): every piece ends at a failure time of its stratum and the
*!                weights are left-continuous, so this is the weight AT the failure time.
*! evalat(mid)    midpoint; evalat(start) piece start. For banded splits (bandpoints() > 1) the weight varies inside a piece, any
*!                evaluation point is an approximation: compare end / mid / start as a sensitivity analysis.
*!
*! TRUNCATION  truncate(lo hi)   percentiles, e.g. truncate(1 99); creates <generate>_pswt = min(max(psw, q_lo), q_hi) and makes it r(weightvar).
*!   The percentiles are PERSON-TIME WEIGHTED over the pieces used in the analysis (adherent pieces if adherence() is given, otherwise all pieces with a
*!   weight): each piece counts with its length time1-time0, so the bounds do not depend on how finely the rows were split (m = 1 or bands).
*!   Definition: q_p = smallest weight w such that the share of person-time with weight <= w is at least p/100 (lower quantile, no interpolation).
*! DIAGNOSTICS  (person-time weighted, adherent pieces if adherence() is given): mean, min, p1, p50, p99, max, effective sample size as a fraction
*!   of person-time  ESS = (sum w t)^2 / (sum t * sum w^2 t);  overall and, with arm(varname), separately for each level of arm.
*!
*! WEIGHT WARNING (0.3.0; thresholds from the simulations of Session 5a/5b: ESS 40-55% = stable (S3, S9), 19% = acceptable after truncation (S12), 6% and 1% = unstable
*!   untruncated estimates with SD 3x the unweighted SD and under-covering SEs (S10, S11)).  The smallest ESS over the overall and the by-arm values (untruncated weight) is compared with
*!     esscaution(#)  percent, default 20:  below it a CAUTION is printed (r(warn)=1)
*!     esssevere(#)   percent, default 5:   below it a SEVERE warning is printed (r(warn)=2)
*!     maxwarn(#)     default 100:  largest weight above it gives a note (r(warn_max)=1)
*!   nowarn suppresses the printed text only (r(warn), r(warn_max) are still returned).  Thresholds are heuristics, not tests: they cannot detect unmeasured confounding.
*!   With truncate() the ESS of the truncated weight is also shown and returned (r(ess_trunc)).
*!
*! Weights are positive for usable rows; zero-length rows have missing parts and therefore missing weights (stset ignores them).
*! The log weight is linear in t inside a row, so  log W(mid) = (log W(start) + log W(end)) / 2  and consecutive pieces join continuously.
*!
*! SYNTAX
*!   ipwadh_eval [if] [in], TIME0(varname) TIME1(varname) [ EVALat(end|mid|start) PREfix(name) GENerate(name) ADHerence(varname)
*!                           ARM(varname) TRUNCate(# #) ESSCaution(#) ESSSevere(#) MAXWarn(#) NOWarn REPLACE ]
*!   prefix(): prefix of the variables made by ipwadhc/ipwadh_split (default _ipw): <prefix>_orig <prefix>_piece <prefix>_E0 _J0 _dA _Ec0 _dC
*!   generate(): prefix of the new variables (default = prefix): <generate>_pswa (adherence), <generate>_pswc (death), <generate>_psw (product), <generate>_pswt (truncated)
*!   time0()/time1() are the PIECE start and end (the variables ipwadh_split has already cut); the data order is restored afterwards.
*!   Returns r(N), r(n_missing), r(evalat), r(weightvar), r(mean), r(min), r(p1), r(p50), r(p99), r(max), r(ess), r(pt), r(n_used), r(ess_min), r(warn), r(warn_max);
*!   with arm(): r(n_levels), r(level_#), r(mean_#), r(p1_#), r(p50_#), r(p99_#), r(max_#), r(ess_#), r(pt_#), r(n_#)  (# = 1, 2, ... in order of the arm values);
*!   with truncate(): r(trunc_lo), r(trunc_hi), r(share_lo), r(share_hi) (person-time shares set to the bounds), r(ess_trunc).

program define ipwadh_eval, rclass sortpreserve
    version 15
    syntax [if] [in], TIME0(varname) TIME1(varname) [ EVALat(string) PREfix(name) GENerate(name) ADHerence(varname) ///
        ARM(varname) TRUNCate(numlist min=2 max=2 >=0 <=100) ESSCaution(real 20) ESSSevere(real 5) MAXWarn(real 100) NOWarn DETail REPLACE ]

    if "`prefix'" == "" local prefix _ipw
    if "`generate'" == "" local generate `prefix'
    if "`evalat'" == "" local evalat end
    if !inlist("`evalat'", "end", "mid", "start") {
        di as err "evalat() must be end, mid or start"
        exit 198
    }
    if `esssevere' <= 0 | `esscaution' > 100 | `esssevere' >= `esscaution' {
        di as err "esssevere(#) < esscaution(#) required, both in (0, 100] percent"
        exit 198
    }
    if `maxwarn' <= 0 {
        di as err "maxwarn(#) must be positive"
        exit 198
    }
    if "`truncate'" != "" {
        local plo : word 1 of `truncate'
        local phi : word 2 of `truncate'
        if `plo' >= `phi' {
            di as err "truncate(lo hi): lo must be smaller than hi"
            exit 198
        }
    }
    foreach s in orig piece E0 J0 dA Ec0 dC {
        capture confirm numeric variable `prefix'_`s'
        if _rc {
            di as err "variable `prefix'_`s' not found (or not numeric)."
            di as err "Run  ipwadhc, ... parts  on the unsplit episodes, then  ipwadh_split, ... , and only then ipwadh_eval."
            exit 111
        }
    }
    local outs pswa pswc psw
    if "`truncate'" != "" local outs `outs' pswt
    foreach s of local outs {
        capture confirm variable `generate'_`s'
        if !_rc {
            if "`replace'" == "" {
                di as err "variable `generate'_`s' already exists; use option replace"
                exit 110
            }
            drop `generate'_`s'
        }
    }

    marksample touse, novarlist
    markout `touse' `time0' `time1' `prefix'_orig `prefix'_piece
    quietly count if `touse'
    if r(N) == 0 error 2000
    local nsel = r(N)

    * pieces of one original row must all be selected, adjacent after sorting, numbered 1..k and contiguous in time
    tempvar ok
    sort `prefix'_orig `prefix'_piece, stable
    quietly by `prefix'_orig: gen byte `ok' = (`prefix'_piece == _n)
    quietly count if `touse' & !`ok'
    if r(N) > 0 {
        di as err "`prefix'_piece is not 1,2,3,... within `prefix'_orig on " r(N) " rows: the data are not the output of ipwadh_split (rows added, dropped or changed?)"
        exit 459
    }
    tempvar gap
    quietly by `prefix'_orig: gen byte `gap' = (_n > 1) & abs(`time0' - `time1'[_n-1]) > 1e-9 * max(1, abs(`time0'))
    quietly count if `touse' & `gap'
    if r(N) > 0 {
        di as err "pieces are not contiguous in time within `prefix'_orig on " r(N) " rows (time0() must equal the previous piece's time1())"
        exit 459
    }
    quietly count if `touse' & `time1' < `time0'
    if r(N) > 0 {
        di as err "time1() < time0() on " r(N) " rows"
        exit 459
    }

    * time since the start of the original row
    tempvar s tt
    quietly by `prefix'_orig: gen double `s' = `time0'[1]
    if "`evalat'" == "end"        quietly gen double `tt' = `time1' - `s'
    else if "`evalat'" == "start" quietly gen double `tt' = `time0' - `s'
    else                          quietly gen double `tt' = (`time0' + `time1') / 2 - `s'

    quietly gen double `generate'_pswa = exp(`prefix'_E0 - `prefix'_dA * `tt' + `prefix'_J0)  if `touse'
    quietly gen double `generate'_pswc = exp(`prefix'_Ec0 - `prefix'_dC * `tt')               if `touse'
    quietly gen double `generate'_psw  = `generate'_pswa * `generate'_pswc                    if `touse'
    label variable `generate'_pswa "IPW adherence at piece `evalat' (continuous time)"
    label variable `generate'_pswc "IPW death/censoring at piece `evalat' (continuous time)"
    label variable `generate'_psw  "Combined IPW at piece `evalat' (untruncated)"

    quietly count if `touse' & missing(`generate'_psw)
    local nmiss = r(N)

    * ---- the pieces used in the analysis (dg), their lengths, person-time weighted diagnostics ---------------------------------
    tempvar dg len
    if "`adherence'" != "" quietly gen byte `dg' = `touse' & `adherence' == 1 & !missing(`generate'_psw) & `time1' > `time0'
    else                   quietly gen byte `dg' = `touse' & !missing(`generate'_psw) & `time1' > `time0'
    quietly gen double `len' = `time1' - `time0'
    quietly count if `dg'
    if r(N) == 0 {
        di as err "no pieces with a weight and positive length"
        exit 2000
    }
    local wv `generate'_psw
    local ess_t = .

    * truncation: person-time weighted percentiles of the combined weight over dg, applied to every selected piece
    if "`truncate'" != "" {
        _ipwadh_wq `generate'_psw `len' `dg' "`plo' `phi'"
        tempname qlo qhi                   // scalars, not locals: a local keeps ~16 digits, the bounds must be the exact weights
        scalar `qlo' = r(q1)
        scalar `qhi' = r(q2)
        quietly gen double `generate'_pswt = min(max(`generate'_psw, `qlo'), `qhi') if `touse' & !missing(`generate'_psw)
        label variable `generate'_pswt "Combined IPW at piece `evalat', truncated at person-time p`plo'/p`phi'"
        quietly summarize `len' if `dg', meanonly
        local pt_all = r(sum)
        quietly summarize `len' if `dg' & `generate'_psw < `qlo', meanonly
        local share_lo = cond(r(N) > 0, r(sum), 0) / `pt_all'
        quietly summarize `len' if `dg' & `generate'_psw > `qhi', meanonly
        local share_hi = cond(r(N) > 0, r(sum), 0) / `pt_all'
        * effective sample size of the TRUNCATED weight
        tempvar tl tl2
        quietly gen double `tl'  = `generate'_pswt * `len' if `dg'
        quietly gen double `tl2' = `generate'_pswt^2 * `len' if `dg'
        quietly summarize `tl', meanonly
        local s1_t = r(sum)
        quietly summarize `tl2', meanonly
        local ess_t = `s1_t'^2 / (`pt_all' * r(sum))
        local wv `generate'_pswt
    }

    if "`detail'" != "" {
        di
        di as txt "ipwadh_eval: weights attached to the pieces, evaluated at the piece " as res "`evalat'"
        di as txt "{hline 78}"
        di as txt "  pieces selected / without weight (zero-length rows)" _col(56) as res %8.0fc `nsel' as txt " / " as res %6.0fc `nmiss'
        di as txt "  diagnostics of " as res "`generate'_psw" as txt " (untruncated), PERSON-TIME weighted" cond("`adherence'" != "", ", adherent pieces", "")
        di as txt "  " _col(14) "pieces" _col(24) "person-time" _col(38) "mean" _col(46) "p1" _col(54) "p50" _col(62) "p99" _col(70) "max" _col(77) "ESS%"
    }
    _ipwadh_diag `generate'_psw `len' `dg' "all" "`detail'"
    foreach x in n pt mean min p1 p50 p99 max ess {
        local a_`x' = r(`x')
    }

    if "`arm'" != "" {
        quietly levelsof `arm' if `dg', local(levs)
        local nl : word count `levs'
        local i 0
        foreach l of local levs {
            local ++i
            tempvar da
            quietly gen byte `da' = `dg' & `arm' == `l'
            _ipwadh_diag `generate'_psw `len' `da' "`arm'=`l'" "`detail'"
            foreach x in n pt mean p1 p50 p99 max ess {
                local b`i'_`x' = r(`x')
            }
            local lev`i' = `l'
            drop `da'
        }
    }

    * ---- weight warning: smallest ESS over the overall and by-arm values ----
    local essmin = `a_ess'
    if "`arm'" != "" {
        forvalues i = 1/`nl' {
            if `b`i'_ess' < `essmin' local essmin = `b`i'_ess'
        }
    }
    local warn = 0
    if 100 * `essmin' < `esssevere'       local warn = 2
    else if 100 * `essmin' < `esscaution' local warn = 1
    local warnmax = (`a_max' > `maxwarn')

    * ---- returned results (copied from locals; no r() is read after this point) ----
    return scalar n_used = `a_n'
    return scalar pt = `a_pt'
    return scalar mean = `a_mean'
    return scalar min = `a_min'
    return scalar p1 = `a_p1'
    return scalar p50 = `a_p50'
    return scalar p99 = `a_p99'
    return scalar max = `a_max'
    return scalar ess = `a_ess'
    return scalar ess_min = `essmin'
    return scalar warn = `warn'
    return scalar warn_max = `warnmax'
    if "`arm'" != "" {
        return scalar n_levels = `nl'
        forvalues i = 1/`nl' {
            return scalar level_`i' = `lev`i''
            foreach x in n pt mean p1 p50 p99 max ess {
                return scalar `x'_`i' = `b`i'_`x''
            }
        }
    }
    if "`detail'" != "" di as txt "  ESS% = (sum w t)^2 / (sum t * sum w^2 t): effective share of the person-time (100% = all weights equal)"
    if "`truncate'" != "" {
        if "`detail'" != "" {
            di as txt "  truncation at person-time p" as res "`plo'" as txt "/p" as res "`phi'" as txt ":  bounds " as res %9.5f `qlo' as txt " and " as res %9.5f `qhi' ///
               as txt ";  person-time set to the bounds: " as res %6.3f 100*`share_lo' as txt "% below, " as res %6.3f 100*`share_hi' as txt "% above"
            di as txt "  ESS of the truncated weight: " as res %5.1f 100*`ess_t' as txt "%"
            di as txt "  truncated weight: " as res "`generate'_pswt"
        }
        return scalar trunc_lo = `qlo'
        return scalar trunc_hi = `qhi'
        return scalar share_lo = `share_lo'
        return scalar share_hi = `share_hi'
        return scalar ess_trunc = `ess_t'
    }
    if "`detail'" != "" di as txt "{hline 78}"
    else {
        di
                di as txt "ipwadh_eval: " as res %-9.0fc `nsel' as txt " pieces" _c
        if `nmiss' > 0 di as txt " (" as res `nmiss' as txt " zero-length without weight)" _c
        di as txt ";  weight " as res "`generate'_psw" as txt " at piece " as res "`evalat'" as txt cond("`adherence'" != "", ", adherent pieces: ", ", all pieces: ") as res %-9.0fc `a_n' as txt " (" as res %9.0fc `a_pt' as txt " person-time)"
        di as txt "  weight, person-time weighted:  mean " as res %5.3f `a_mean' as txt "  p1 " as res %6.4f `a_p1' as txt "  p50 " as res %5.3f `a_p50' as txt "  p99 " as res %5.2f `a_p99' as txt "  max " as res %6.2f `a_max'
        di as txt "  effective sample size " as res %4.1f 100*`a_ess' as txt "%" _c
        if "`arm'" != "" di as txt " (smallest by arm " as res %4.1f 100*`essmin' as txt "%)" _c
        di as txt ""
        if "`truncate'" != "" {
            di as txt "  truncated at p" as res "`plo'" as txt "/p" as res "`phi'" as txt " (bounds " as res %6.4f `qlo' as txt ", " as res %5.2f `qhi' as txt "; " as res %4.2f 100*(`share_lo' + `share_hi') ///
               as txt "% of the person-time at the bounds):  " as res "`generate'_pswt" as txt ", ESS " as res %4.1f 100*`ess_t' as txt "%"
        }
    }
    if "`nowarn'" == "" {
        if `warn' == 2 {
            di as err "WARNING (severe): effective sample size " %4.1f 100*`essmin' "% (< `esssevere'%): the weighted estimate is unstable and its SE too small."
            di as err "  Simplify the adherence model (fewer / coarser predictors, look for near-violations of positivity) and report the unweighted, untruncated and truncated analyses side by side."
        }
        else if `warn' == 1 {
            di as txt "CAUTION: effective sample size " as res %4.1f 100*`essmin' as txt "% (< `esscaution'%): weights are heavy-tailed. Compare the unweighted, untruncated and truncated (truncate(1 99)) analyses."
        }
        if `warnmax' {
            di as txt "NOTE: largest weight " as res %9.2f `a_max' as txt " (> `maxwarn'): check the pieces with the largest weights."
        }
        if "`detail'" != "" & `warn' > 0 & "`truncate'" != "" {
            di as txt "  After truncation the ESS is " as res %4.1f 100*`ess_t' as txt "%, but truncation also removes part of the correction (bias towards the unweighted estimate)."
        }
    }
    if "`evalat'" != "end" di as txt "NOTE: evalat(`evalat') is an approximation unless pieces are very short; evalat(end) is exact for the default split."
    if "`detail'" != "" {
        di as txt "Analyse with:  stset `time1' [pw=`wv'] if <adherent>, time0(`time0') failure()    (no id(): time-varying weights)"
        di as txt "              stcox ..., breslow vce(cluster id)      (pweights give robust SEs; iweights are not allowed with robust)"
    }

    return scalar N = `nsel'
    return scalar n_missing = `nmiss'
    return local evalat "`evalat'"
    return local weightvar "`wv'"
end


* person-time weighted lower quantiles of key over the observations with touse == 1; p-list in percent; returns r(q1), r(q2), ...
* (sorts the data; the caller runs under sortpreserve)
program define _ipwadh_wq, rclass
    args key len touse plist
    tempvar k l c
    quietly gen double `k' = `key' if `touse'
    quietly gen double `l' = `len' if `touse'
    sort `k'
    quietly gen double `c' = sum(`l')
    quietly summarize `l', meanonly
    local tot = r(sum)
    local i 0
    foreach p of local plist {
        local ++i
        quietly summarize `k' if !missing(`k') & `c' >= `tot' * `p' / 100 * (1 - 1e-12), meanonly
        return scalar q`i' = r(min)
    }
end


* one printed diagnostics row (person-time weighted); returns r(n) r(pt) r(mean) r(min) r(p1) r(p50) r(p99) r(max) r(ess)
program define _ipwadh_diag, rclass
    args w len touse label show
    tempvar wl w2l
    quietly gen double `wl'  = `w' * `len'    if `touse'
    quietly gen double `w2l' = `w'^2 * `len'  if `touse'
    quietly summarize `len' if `touse', meanonly
    local n = r(N)
    local pt = r(sum)
    quietly summarize `wl' if `touse', meanonly
    local sw = r(sum)
    quietly summarize `w2l' if `touse', meanonly
    local sw2 = r(sum)
    quietly summarize `w' if `touse', meanonly
    local mn = r(min)
    local mx = r(max)
    local mean = `sw' / `pt'
    local ess = `sw'^2 / (`pt' * `sw2')
    _ipwadh_wq `w' `len' `touse' "1 50 99"
    local p1 = r(q1)
    local p50 = r(q2)
    local p99 = r(q3)
    if "`show'" != "" di as txt "  " abbrev("`label'", 11) _col(13) as res %7.0fc `n' _col(22) %12.2f `pt' _col(36) %7.4f `mean' _col(44) %7.4f `p1' _col(52) %7.4f `p50' ///
       _col(60) %7.4f `p99' _col(68) %7.3f `mx' _col(77) %4.1f 100*`ess'
    return scalar n = `n'
    return scalar pt = `pt'
    return scalar mean = `mean'
    return scalar min = `mn'
    return scalar p1 = `p1'
    return scalar p50 = `p50'
    return scalar p99 = `p99'
    return scalar max = `mx'
    return scalar ess = `ess'
end
