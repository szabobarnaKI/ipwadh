* ipwadh_mata.do -- Mata support functions for ipwadh_split.ado (loaded on demand; lmbuild packaging later)
* NOT YET RUN IN STATA (written in a container without Stata) -- see LOG.md
* Method: exact-minimum, budgeted thinning.  See ipwadh_handoff_plan.md section 1c.
* v0.2.2: drops any earlier definitions first, so the file can be run repeatedly in one Stata session
*         (fixes "ipwadh_cnt() already exists", r(3000)).
version 15
capture mata: mata drop ipwadh_cnt()
capture mata: mata drop ipwadh_cutpoints()
capture mata: mata drop ipwadh_segs()
capture mata: mata drop ipwadh_kvec()
capture mata: mata drop ipwadh_rows()
capture mata: mata drop ipwadh_auto()
capture mata: mata drop ipwadh_prep()
capture mata: mata drop ipwadh_fill()
mata:
mata set matastrict on

// number of elements of the sorted segment B[first .. first+n-1] that are < x (strict=1) or <= x (strict=0)
real scalar ipwadh_cnt(real colvector B, real scalar first, real scalar n, real scalar x, real scalar strict)
{
    real scalar lo, hi, mid, v
    lo = 0
    hi = n
    while (lo < hi) {
        mid = floor((lo + hi + 1) / 2)
        v = B[first + mid - 1]
        if ((strict & v < x) | (!strict & v <= x)) lo = mid
        else hi = mid - 1
    }
    return(lo)
}

// Cut points as a matrix BM = (stratum, time), sorted by stratum then time.
//  failure points FP = distinct (stratum, failure time) pairs (adherent rows only if adhonly)
//  C = distinct union of all failure times and all adherence-change times (chgv; missing = no change)
//  m == 1 : exact method, boundaries = FP (cutall: C for every stratum)
//  m  > 1 : groups of m consecutive points of C; each stratum gets, per group, ONE boundary at its LAST OWN failure
//           time in the group (none if it has no failure there) => boundaries are a subset of FP, rows(m) <= rows(1)
// sets locals nC (=|C|), nFp (=rows(FP)), nB (=rows(BM))
real matrix ipwadh_cutpoints(string scalar t1v, string scalar fv, string scalar advar, string scalar chgv,
                             string scalar sgv, string scalar touse, real scalar adhonly, real scalar m,
                             real scalar cutall, real scalar S)
{
    real colvector t1, f, sg, sel, A, C, Cc, j, g, lastf
    real matrix FP, BM
    real scalar i, n, s
    t1 = st_data(., t1v, touse)
    f  = st_data(., fv,  touse)
    sg = st_data(., sgv, touse)
    sel = (f :== 1)
    if (adhonly) sel = sel :& (st_data(., advar, touse) :== 1)
    FP = select((sg, t1), sel)
    if (rows(FP) > 0) FP = uniqrows(sort(FP, (1, 2)))
    A = st_data(., chgv, touse)
    A = select(A, A :< .)
    if (rows(FP) > 0) Cc = FP[., 2] \ A
    else Cc = A
    st_local("nFp", strofreal(rows(FP), "%18.0f"))
    if (rows(Cc) == 0) {
        st_local("nC", "0")
        st_local("nB", "0")
        return(J(0, 2, .))
    }
    C = uniqrows(sort(Cc, 1))
    st_local("nC", strofreal(rows(C), "%18.0f"))
    if (m == 1) {
        if (cutall) {
            BM = J(0, 2, .)
            for (s = 1; s <= S; s++) BM = BM \ (J(rows(C), 1, s), C)
        }
        else BM = FP
    }
    else if (rows(FP) == 0) BM = J(0, 2, .)
    else {
        n = rows(FP)
        j = J(n, 1, .)
        for (i = 1; i <= n; i++) j[i] = ipwadh_cnt(C, 1, rows(C), FP[i, 2], 0)   // 1-based position of the time in C
        g = floor((j :- 1) :/ m)
        lastf = J(n, 1, 1)
        for (i = 1; i < n; i++) lastf[i] = (FP[i, 1] != FP[i+1, 1]) | (g[i] != g[i+1])
        BM = select(FP, lastf)
    }
    st_local("nB", strofreal(rows(BM), "%18.0f"))
    return(BM)
}

// first position and count of each stratum's segment in BM
void ipwadh_segs(real matrix BM, real scalar S, real colvector sfirst, real colvector scount)
{
    real scalar i, s
    sfirst = J(S, 1, 0)
    scount = J(S, 1, 0)
    for (i = 1; i <= rows(BM); i++) {
        s = BM[i, 1]
        if (scount[s] == 0) sfirst[s] = i
        scount[s] = scount[s] + 1
    }
}

// per row: lo = #cut points of its stratum <= t0 ; k = # cut points strictly inside (t0,t1). Returns rows after the split.
real scalar ipwadh_kvec(real matrix BM, real colvector t0, real colvector t1, real colvector a, real colvector sg,
                        real scalar adhonly, real scalar S, real colvector k, real colvector lo)
{
    real colvector sfirst, scount, Bt
    real scalar i, n, s, h
    n  = rows(t0)
    k  = J(n, 1, 0)
    lo = J(n, 1, 0)
    ipwadh_segs(BM, S, sfirst, scount)
    if (rows(BM) > 0) {
        Bt = BM[., 2]
        for (i = 1; i <= n; i++) {
            if (adhonly & a[i] != 1) continue
            s = sg[i]
            if (scount[s] == 0) continue
            lo[i] = ipwadh_cnt(Bt, sfirst[s], scount[s], t0[i], 0)
            h     = ipwadh_cnt(Bt, sfirst[s], scount[s], t1[i], 1)
            k[i]  = max((h - lo[i], 0))
        }
    }
    return(n + sum(k))
}

// rows after the split of the selected rows for a given m
real scalar ipwadh_rows(real colvector t0, real colvector t1, real colvector a, real colvector sg, string scalar t1v,
                        string scalar fv, string scalar advar, string scalar chgv, string scalar sgv, string scalar touse,
                        real scalar adhonly, real scalar m, real scalar cutall, real scalar S)
{
    real matrix BM
    real colvector k, lo
    BM = ipwadh_cutpoints(t1v, fv, advar, chgv, sgv, touse, adhonly, m, cutall, S)
    return(ipwadh_kvec(BM, t0, t1, a, sg, adhonly, S, k, lo))
}

// Choose m. Always reports rows_exact (m=1), nC, mmax.  budget = max rows (selected rows after split) allowed.
// m_sel = 1 if the exact split fits; else smallest m (doubling, then bisection) that fits; m_sel = 0 if even mmax does not fit.
void ipwadh_auto(string scalar t0v, string scalar t1v, string scalar fv, string scalar advar, string scalar chgv,
                 string scalar sgv, string scalar touse, real scalar adhonly, real scalar cutall, real scalar S,
                 real scalar budget)
{
    real colvector t0, t1, a, sg
    real scalar n1, mmax, lo, hi, mid, rows_hi, nC
    t0 = st_data(., t0v, touse)
    t1 = st_data(., t1v, touse)
    a  = st_data(., advar, touse)
    sg = st_data(., sgv, touse)
    n1 = ipwadh_rows(t0, t1, a, sg, t1v, fv, advar, chgv, sgv, touse, adhonly, 1, cutall, S)
    nC = strtoreal(st_local("nC"))
    mmax = max((1, floor(nC / 2)))
    st_local("rows_exact", strofreal(n1, "%18.0f"))
    st_local("mmax", strofreal(mmax, "%18.0f"))
    if (n1 <= budget) {
        st_local("m_sel", "1")
        return
    }
    lo = 1
    hi = 2
    while (hi < mmax) {
        if (ipwadh_rows(t0, t1, a, sg, t1v, fv, advar, chgv, sgv, touse, adhonly, hi, cutall, S) <= budget) break
        lo = hi
        hi = 2 * hi
    }
    if (hi > mmax) hi = mmax
    rows_hi = ipwadh_rows(t0, t1, a, sg, t1v, fv, advar, chgv, sgv, touse, adhonly, hi, cutall, S)
    if (rows_hi > budget) {
        st_local("m_sel", "0")
        st_local("rows_floor", strofreal(rows_hi, "%18.0f"))
        return
    }
    while (hi - lo > 1) {
        mid = floor((lo + hi) / 2)
        if (ipwadh_rows(t0, t1, a, sg, t1v, fv, advar, chgv, sgv, touse, adhonly, mid, cutall, S) <= budget) hi = mid
        else lo = mid
    }
    st_local("m_sel", strofreal(hi, "%18.0f"))
}

// fill per-row k and lo for a chosen m
void ipwadh_prep(string scalar t0v, string scalar t1v, string scalar fv, string scalar advar, string scalar chgv,
                 string scalar sgv, string scalar touse, string scalar kv, string scalar lov, real scalar adhonly,
                 real scalar m, real scalar cutall, real scalar S)
{
    real matrix BM
    real colvector t0, t1, a, sg, k, lo
    BM = ipwadh_cutpoints(t1v, fv, advar, chgv, sgv, touse, adhonly, m, cutall, S)
    t0 = st_data(., t0v, touse)
    t1 = st_data(., t1v, touse)
    a  = st_data(., advar, touse)
    sg = st_data(., sgv, touse)
    (void) ipwadh_kvec(BM, t0, t1, a, sg, adhonly, S, k, lo)
    st_store(., kv,  touse, k)
    st_store(., lov, touse, lo)
}

// after -expand-: piece j (1..k+1) of a row gets its own (t0,t1) cut at its stratum's boundaries
void ipwadh_fill(string scalar t0v, string scalar t1v, string scalar fv, string scalar advar, string scalar chgv,
                 string scalar sgv, string scalar touse, string scalar piecev, string scalar kv, string scalar lov,
                 real scalar adhonly, real scalar m, real scalar cutall, real scalar S)
{
    real matrix BM
    real colvector t0, t1, j, k, lo, sg, n0, n1, sfirst, scount, Bt
    real scalar i, n, s, f0
    BM = ipwadh_cutpoints(t1v, fv, advar, chgv, sgv, touse, adhonly, m, cutall, S)  // copies carry identical t1/failure/chgv
    t0 = st_data(., t0v, touse)
    t1 = st_data(., t1v, touse)
    j  = st_data(., piecev, touse)
    k  = st_data(., kv, touse)
    lo = st_data(., lov, touse)
    sg = st_data(., sgv, touse)
    ipwadh_segs(BM, S, sfirst, scount)
    n0 = t0
    n1 = t1
    n  = rows(t0)
    if (rows(BM) > 0) {
        Bt = BM[., 2]
        for (i = 1; i <= n; i++) {
            if (k[i] > 0) {
                s  = sg[i]
                f0 = sfirst[s]
                if (j[i] > 1)     n0[i] = Bt[f0 + lo[i] + j[i] - 2]
                if (j[i] <= k[i]) n1[i] = Bt[f0 + lo[i] + j[i] - 1]
            }
        }
    }
    st_store(., t0v, touse, n0)
    st_store(., t1v, touse, n1)
}
end
