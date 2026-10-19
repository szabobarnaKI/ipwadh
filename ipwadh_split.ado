*! ipwadh_split 0.2.2  Exact-minimum split of episode data at failure times, with optional budgeted thinning (bands).
*! 0.2.2: compact default output (one line after the split); the size table is shown with option detail or check, for banded splits and whenever a memory warning applies. Numbers unchanged.
*! 0.2.1 passed test_ipwadh_split.do in Stata 19.
*!
*! WHAT IT DOES
*!   Cuts rows at failure times lying strictly inside (time0,time1), like -stsplit, at(failures)-, but only where an
*!   on-treatment analysis (-stset ... [iw=w] if adherence==1-, optionally stratified) NEEDS a cut:
*!     - only ADHERENT rows are cut (allrows: all rows, as stsplit does);
*!     - only at failure times of adherent rows;
*!     - with strata(varlist): a row is cut only at failure times of ITS OWN stratum (strata must be constant within a row);
*!     - never at other persons' adherence changes (rows already start/end at each person's own changes).
*!   This default, bandpoints(1), is exact (identical Breslow estimates and robust SEs, verified in Python on S3-S7, also stratified).
*!
*! BANDS (approximation; only when the exact split is too big)
*!   bandpoints(#):  # = 1 (default, exact) or an integer 2..50% of |C|. C = union of failure times and adherence-change times.
*!                   Groups of # consecutive points of C; each stratum gets ONE boundary per group, at its LAST OWN FAILURE time
*!                   in the group (none if no failure). Boundaries are a subset of the exact cuts, so rows never exceed the exact count.
*!   bandpoints(auto): m = 1 if the exact split fits the row budget, otherwise the smallest m that fits (doubling + bisection;
*!                   counting passes only, data untouched). Budget = min(maxrows(), maxgb()-based, Stata edition/max_memory limits).
*!   Weights are constant within a piece: evaluate them at the piece midpoint (most stable in tests), start or end.
*!
*! SYNTAX
*!   ipwadh_split [if] [in], id(varname) time0(varname) time1(varname) adherence(varname) failure(varname)
*!        [ strata(varlist) allrows lastonly(varlist) bandpoints(# | auto) cutall check detail maxgb(#) maxrows(#) force maxobs(#) generate(name) ]
*!   lastonly(): indicator variables that must be 1 only on the last piece of a row (failure is always included; add dead/ltfu).
*!   check: report only, do not split (prints the size table).  detail: print the size table also for the split itself.  maxgb(#): memory limit (default 8 GB); maxrows(#): row limit.
*!   force: proceed with the split although the ESTIMATED MEMORY need exceeds maxgb(), maxrows() or Stata's -set max_memory-
*!          (to test whether Stata/Mata and the operating system, e.g. Windows paging/virtual memory, can cope; with bandpoints(auto)
*!          force removes the memory-based budget, so auto then returns m=1 unless the observation limit is hit).
*!          force NEVER overrides the maximum-observations limit of the Stata edition (r(900), Stata could not hold the rows).
*!          If -set max_memory- is finite, Stata itself may still refuse the allocation: -set max_memory .- first.
*!          The split is attempted; if Stata's expand fails (r(900)/r(901)/r(909)/r(950)) the command reports it.
*!   maxobs(#): testing aid, pretends the edition allows at most # observations (to test the hard limit).
*!   cutall: with bandpoints(1) keep cuts at every point of C (redundant, many more rows; for testing).
*!   Creates <generate>_orig (original observation number) and <generate>_piece (1..k+1); default prefix _ipw.
*!   Data are left in the original order, pieces adjacent. Rows outside if/in are untouched.
*!   Returns r(exact) (1 if m==1), r(m), r(N_exact), r(N_after), r(n_cutpoints), r(max_bandpoints), r(est_gb), r(ok).

program define ipwadh_split, rclass
    version 15
    syntax [if] [in], ID(varname) TIME0(varname) TIME1(varname) ADHerence(varname) FAILure(varname) ///
        [ STRata(varlist) ALLRows LASTonly(varlist) BANDpoints(string) CUTALL CHECK DETail ///
          MAXGB(real 8) MAXRows(real 0) FORCE MAXOBS(real 0) GENerate(name) ]

    capture mata: mata which ipwadh_auto()
    if _rc {
        quietly findfile ipwadh_mata.do
        quietly run `"`r(fn)'"'
    }
    if "`generate'" == "" local generate _ipw
    foreach s in orig piece {
        capture confirm variable `generate'_`s'
        if !_rc {
            di as err "variable `generate'_`s' already exists"
            exit 110
        }
    }
    local adhonly = ("`allrows'" == "")
    local cutall = ("`cutall'" != "")
    local auto 0
    local m 1
    if "`bandpoints'" == "auto" local auto 1
    else if "`bandpoints'" != "" {
        capture confirm integer number `bandpoints'
        if _rc | real("`bandpoints'") < 1 {
            di as err "bandpoints() must be a positive integer or auto"
            exit 198
        }
        local m = real("`bandpoints'")
    }

    marksample touse, novarlist
    markout `touse' `id' `time0' `time1' `adherence' `failure'
    if "`strata'" != "" markout `touse' `strata', strok
    quietly count if `touse' & !inlist(`adherence', 0, 1)
    if r(N) > 0 {
        di as err "adherence() must be coded 0/1"
        exit 198
    }
    quietly count if `touse' & !inlist(`failure', 0, 1)
    if r(N) > 0 {
        di as err "failure() must be coded 0/1"
        exit 198
    }
    quietly count if `touse' & `time1' < `time0'
    if r(N) > 0 {
        di as err "time1() < time0() on " r(N) " rows"
        exit 459
    }
    quietly count if `touse'
    local nsel = r(N)
    if `nsel' == 0 error 2000

    * ---------- strata index, adherence-change times (data order restored afterwards) ----------
    tempvar sg ord chgt kk lo
    if "`strata'" == "" quietly gen long `sg' = 1 if `touse'
    else quietly egen long `sg' = group(`strata') if `touse'
    quietly summarize `sg' if `touse', meanonly
    local S = r(max)
    quietly gen long `ord' = _n
    quietly gen double `chgt' = .
    sort `touse' `id' `time0' `time1' `ord'
    quietly by `touse' `id': replace `chgt' = `time1' if `touse' & _n < _N & `adherence' != `adherence'[_n+1]
    sort `ord'

    * ---------- limits of this Stata, memory budget ----------
    local N0 = _N
    quietly describe, short
    local width = r(width)
    local maxN = .
    capture local maxN = c(max_N_theory)
    if missing(`maxN') {
        if c(MP) == 1 local maxN = 1099511627775
        else local maxN = 2147483619
    }
    local edition "unknown"
    capture local edition = c(edition)
    if `maxobs' > 0 local maxN = min(`maxN', `maxobs')     // testing aid
    local maxmem = .
    capture local maxmem = c(max_memory)
    local force_on = ("`force'" != "")
    * peak need per row (heuristic): during the split = data + orig/piece (8) + temporary variables (56) + Mata vectors (64);
    * downstream (stset + weights) needs about width + 66, so the split peak governs
    local bytes_per_row = `width' + 8 + 56 + 64
    local bytes_factor = `bytes_per_row' * 1.25      // 1.25 = working-copy overhead (heuristic)
    local soft_rows = `maxgb' * 1024^3 / `bytes_factor'
    if `maxrows' > 0 local soft_rows = min(`soft_rows', `maxrows')
    local hard_rows = `maxN'
    if !missing(`maxmem') local hard_rows = min(`hard_rows', `maxmem' / `bytes_factor')
    * force removes every memory-based limit from the budget; only the observation limit stays
    if `force_on' {
        local soft_rows = `maxN'
        local hard_rows = `maxN'
    }
    local budget_all = min(`soft_rows', `hard_rows')           // rows allowed in the whole dataset after the split
    local budget_sel = `budget_all' - (`N0' - `nsel')          // rows allowed for the selected rows

    * ---------- exact count, choice of m (nothing is modified; counting passes only) ----------
    local nC 0
    local nFp 0
    local nB 0
    local rows_exact 0
    local mmax 1
    local m_sel 1
    local rows_floor 0
    local bud = cond(`auto', `budget_sel', 1e18)
    mata: ipwadh_auto("`time0'", "`time1'", "`failure'", "`adherence'", "`chgt'", "`sg'", "`touse'", `adhonly', `cutall', `S', `bud')
    local N_exact = `N0' - `nsel' + `rows_exact'
    if `auto' {
        if `m_sel' == 0 {
            local fl = `N0' - `nsel' + `rows_floor'
            di as err "bandpoints(auto): the row budget (" %15.0fc `budget_all' ") cannot be met: even bandpoints(`mmax') needs " %15.0fc `fl' " rows."
            di as err "Splitting NOT started. Raise maxgb()/maxrows(), use strata(), analyse a subsample, or run on a larger machine."
            return scalar ok = 0
            exit 950
        }
        local m = `m_sel'
    }
    if `m' > `mmax' {
        di as err "bandpoints(`m') exceeds 50% of the `nC' cut points (union of failure and adherence-change times); maximum is `mmax'"
        exit 198
    }
    quietly gen double `kk' = . if `touse'
    quietly gen double `lo' = . if `touse'
    mata: ipwadh_prep("`time0'", "`time1'", "`failure'", "`adherence'", "`chgt'", "`sg'", "`touse'", "`kk'", "`lo'", `adhonly', `m', `cutall', `S')
    quietly summarize `kk' if `touse', meanonly
    local ksum = r(sum)
    local kmax = r(max)
    local N1 = `N0' + `ksum'
    local est_gb = `N1' * `bytes_factor' / (1024^3)

    local hard_N = (`N1' > `maxN')                                    // never overridable
    local memlim = (!missing(`maxmem') & `N1' * `bytes_factor' > `maxmem')   // set max_memory: overridable with force
    local soft = (`est_gb' > `maxgb' | (`maxrows' > 0 & `N1' > `maxrows'))   // maxgb()/maxrows(): overridable with force

    * the size table is shown on request (detail, check), for approximate bands and whenever a memory limit is hit
    local showtab = ("`detail'" != "" | "`check'" != "" | `m' > 1 | `hard_N' | `memlim' | `soft')

    local ttl "adherent rows only"
    if !`adhonly' local ttl "all rows"
    if "`strata'" != "" local ttl "`ttl', strata-specific cuts (`S' strata)"
    if `showtab' {
        di
        di as txt "ipwadh_split: pre-split size check  (`ttl')"
        di as txt "{hline 70}"
        di as txt "  rows selected / in memory" _col(46) as res %16.0fc `nsel' as txt " / " as res %14.0fc `N0'
        di as txt "  cut points |C| (failures + adherence changes)" _col(46) as res %16.0fc `nC'
        di as txt "  (stratum, failure time) pairs" _col(46) as res %16.0fc `nFp'
        di as txt "  rows after EXACT split (m=1)" _col(46) as res %16.0fc `N_exact'
        if `m' == 1 {
            di as txt "  method" _col(46) as res "exact (bandpoints(1))"
        }
        else {
            local mean_band = `nFp' / max(`nB', 1)
            local autotxt ""
            if `auto' local autotxt " (chosen by auto)"
            di as txt "  method" _col(46) as res "APPROXIMATE bands, bandpoints(`m')" as txt "`autotxt'"
            di as txt "  boundaries used / mean failure times per band" _col(46) as res %16.0fc `nB' as txt " / " as res %9.1f `mean_band'
        }
        di as txt "  rows after split" _col(46) as res %16.0fc `N1' as txt "   (" as res %5.1f 100*`N1'/`N_exact' as txt "% of exact)"
        di as txt "  max pieces for one row" _col(46) as res %16.0fc (`kmax' + 1)
        di as txt "  Stata edition / max observations" _col(46) as res "`edition'" as txt " / " as res %14.0fc `maxN'
        di as txt "  estimated peak memory during split" _col(46) as res %16.1f `est_gb' as txt " GB (heuristic)"
        di as txt "  set max_memory / maxgb()" _col(46) as res cond(missing(`maxmem'), ".", string(`maxmem'/1024^3, "%9.1f") + " GB") as txt " / " as res %6.1f `maxgb' as txt " GB"
        di as txt "{hline 70}"
    }
    if `m' > 1 {
        di as txt "NOTE: banded split is an approximation: evaluate weights at the piece MIDPOINT and compare with a run using"
        di as txt "      a different bandpoints() (e.g. m/2) as a sensitivity analysis."
    }

    return scalar N_before = `N0'
    return scalar N_after = `N1'
    return scalar N_exact = `N_exact'
    return scalar m = `m'
    return scalar exact = (`m' == 1)
    return scalar n_cutpoints = `nC'
    return scalar n_boundaries = `nB'
    return scalar max_bandpoints = `mmax'
    return scalar est_gb = `est_gb'
    return scalar max_N = `maxN'

    local why ""
    local est_txt = string(`est_gb', "%9.1f")
    local maxN_txt = string(`maxN', "%18.0fc")
    if `hard_N' local why "`why' the split needs more observations than this Stata edition can hold (limit `maxN_txt');"
    if `memlim' local why "`why' the estimated memory (`est_txt' GB) exceeds set max_memory;"
    if `soft' local why "`why' the estimated memory (`est_txt' GB) exceeds maxgb(`maxgb') or maxrows();"
    return scalar forced = 0
    if "`why'" != "" {
        if `hard_N' {
            di as err "ERROR:`why'"
            di as err "The maximum-observations limit cannot be overridden (force does not apply). Use bandpoints(auto) / a larger bandpoints(),"
            di as err "strata(), a subsample, or Stata/MP."
            return scalar ok = 0
            if "`check'" == "" exit 900
            exit
        }
        di as err "WARNING:`why'"
        if "`force'" == "" {
            di as err "Splitting NOT started. Options: force (try anyway; the operating system may page to virtual memory), bandpoints(auto) or a"
            di as err "larger bandpoints() (approximate bands), strata(), a subsample of persons, raise maxgb(), more RAM / Stata/MP."
            return scalar ok = 0
            if "`check'" == "" exit 950
            exit
        }
        di as err "force specified: the split will be ATTEMPTED despite the memory warning (exact method kept: bandpoints(`m'))."
        di as err "If set max_memory is finite Stata may still refuse; try -set max_memory .- . Expect heavy paging if RAM is exceeded."
        return scalar forced = 1
    }
    return scalar ok = 1
    if "`check'" != "" exit

    * ---------- split ----------
    tempvar seq ncopy
    quietly gen long `generate'_orig = _n
    quietly replace `kk' = 0 if !`touse'
    quietly gen double `ncopy' = `kk' + 1
    timer clear 97
    timer on 97
    capture quietly expand `ncopy'
    if _rc {
        local rc_expand = _rc
        di as err "expand failed with r(`rc_expand'): Stata could not allocate the " %15.0fc `N1' " rows (data in memory left as before the split call, but"
        di as err "the _orig variable may remain: -drop `generate'_orig- if it exists). Estimated need `est_txt' GB."
        capture drop `generate'_orig
        exit `rc_expand'
    }
    quietly gen double `seq' = _n      // after expand: copies get distinct, increasing numbers
    quietly bysort `generate'_orig (`seq'): gen long `generate'_piece = _n
    mata: ipwadh_fill("`time0'", "`time1'", "`failure'", "`adherence'", "`chgt'", "`sg'", "`touse'", "`generate'_piece", "`kk'", "`lo'", `adhonly', `m', `cutall', `S')
    foreach v of varlist `failure' `lastonly' {
        quietly replace `v' = 0 if `touse' & `generate'_piece <= `kk'
    }
    sort `generate'_orig `generate'_piece
    label variable `generate'_orig "original observation number (before ipwadh_split)"
    label variable `generate'_piece "piece within original row"
    timer off 97
    quietly timer list 97
    local secs = r(t97)
    if `showtab' di as txt "Split done: " as res %12.0fc _N as txt " rows in " as res %9.1f `secs' as txt " s. Linking variables: `generate'_orig, `generate'_piece"
    else {
        di
        di as txt "ipwadh_split: " as res %-9.0fc `N0' as txt " rows -> " as res %-9.0fc _N as txt " pieces (exact split at the failure times of the adherent rows, " as res %-6.0fc `nC' as txt " cut points, " as res %4.1f `secs' as txt " s)"
    }
end
