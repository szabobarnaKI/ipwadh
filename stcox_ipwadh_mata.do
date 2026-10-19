* stcox_ipwadh_mata.do 0.1.0  Mata functions of stcox_ipwadh (run automatically by stcox_ipwadh when needed; NOT YET RUN IN STATA)
* Mirrors ref_se.py (validated in Python: derivative terms vs numerical differences, rel. error ~1e-6).
*
* Model:  U(beta, theta) = sum_p s_p-type weighted Breslow score;  dU/d log w_p = s_p (score-residual contribution of piece p)
*         log w_p = E0 - dA x + J0 + Ec0 - dC x   (x = evaluation time - start of the original row)
*         d log w_p / d theta = closed form in the Poisson linear predictors  ->  B_f = sum_j kappa_fj z_fj'
*         psi_i = U_i + sum_f B_f V_f S_if ;  V = Ainv (sum_i psi_i psi_i') Ainv * G/(G-1)

version 15
capture mata: mata drop ipwadh_cs() ipwadh_cnt_lt() ipwadh_cnt_le() ipwadh_cox_stratum() ipwadh_se_start() ipwadh_se_clear() ipwadh_se_rows() ipwadh_se_addfit() ipwadh_se_parts() ipwadh_se_final()
mata:
mata set matastrict off

// column-wise running sums of a matrix (runningsum() accepts vectors only)
real matrix ipwadh_cs(real matrix X)
{
    real matrix R
    real scalar j
    R = X
    for (j = 1; j <= cols(X); j++) R[., j] = runningsum(X[., j])
    return(R)
}

// number of elements of the ascending vector S that are strictly smaller than x[i]
real colvector ipwadh_cnt_lt(real colvector S, real colvector x)
{
    real colvector o, res
    real scalar i, j, m, n, xi
    m = rows(S)
    n = rows(x)
    res = J(n, 1, 0)
    if (n == 0 | m == 0) return(res)
    o = order(x, 1)
    j = 0
    for (i = 1; i <= n; i++) {
        xi = x[o[i]]
        while (j < m) {
            if (S[j+1] < xi) j++
            else break
        }
        res[o[i]] = j
    }
    return(res)
}

// number of elements of the ascending vector S that are smaller than or equal to x[i]
real colvector ipwadh_cnt_le(real colvector S, real colvector x)
{
    real colvector o, res
    real scalar i, j, m, n, xi
    m = rows(S)
    n = rows(x)
    res = J(n, 1, 0)
    if (n == 0 | m == 0) return(res)
    o = order(x, 1)
    j = 0
    for (i = 1; i <= n; i++) {
        xi = x[o[i]]
        while (j < m) {
            if (S[j+1] <= xi) j++
            else break
        }
        res[o[i]] = j
    }
    return(res)
}

// Weighted Breslow quantities in ONE stratum. Risk set of event time T: t0 < T <= t1. Returns (by reference)
//   s  n x p   score-residual contribution of every piece (event term minus risk-set term, weights included)
//   A  p x p   observed information   U  p x 1   total score
void ipwadh_cox_stratum(real colvector t0, real colvector t1, real colvector ev, real colvector w, real matrix X,
                        real colvector beta, real matrix s, real matrix A, real colvector U)
{
    real scalar n, p, m, i, a, b
    real colvector evi, T, r, o0, o1, i0, i1, c0, c1, S0, jj, D, dL, Lam, ia, ib, v, S2ab
    real matrix cc0, cc1, S1, xbar, Gam, xbe

    n = rows(t0)
    p = cols(X)
    s = J(n, p, 0)
    A = J(p, p, 0)
    U = J(p, 1, 0)
    evi = selectindex(ev :== 1)
    if (rows(evi) == 0) return
    T = uniqrows(t1[evi])
    m = rows(T)
    r = w :* exp(X * beta)
    o0 = order(t0, 1)
    o1 = order(t1, 1)
    i0 = ipwadh_cnt_lt(t0[o0], T)
    i1 = ipwadh_cnt_lt(t1[o1], T)
    c0 = 0 \ runningsum(r[o0])
    c1 = 0 \ runningsum(r[o1])
    S0 = c0[i0 :+ 1] - c1[i1 :+ 1]
    cc0 = J(1, p, 0) \ ipwadh_cs(r[o0] :* X[o0, .])
    cc1 = J(1, p, 0) \ ipwadh_cs(r[o1] :* X[o1, .])
    S1 = cc0[i0 :+ 1, .] - cc1[i1 :+ 1, .]
    xbar = S1 :/ S0
    jj = ipwadh_cnt_lt(T, t1[evi]) :+ 1
    D = J(m, 1, 0)
    for (i = 1; i <= rows(evi); i++) D[jj[i]] = D[jj[i]] + w[evi[i]]
    U = (colsum(w[evi] :* X[evi, .]) - colsum(D :* xbar))'
    for (a = 1; a <= p; a++) {
        for (b = a; b <= p; b++) {
            v = r :* X[., a] :* X[., b]
            c0 = 0 \ runningsum(v[o0])
            c1 = 0 \ runningsum(v[o1])
            S2ab = c0[i0 :+ 1] - c1[i1 :+ 1]
            A[a, b] = sum(D :* (S2ab :/ S0 - xbar[., a] :* xbar[., b]))
            A[b, a] = A[a, b]
        }
    }
    dL = D :/ S0
    Lam = 0 \ runningsum(dL)
    Gam = J(1, p, 0) \ ipwadh_cs(dL :* xbar)
    ia = ipwadh_cnt_le(T, t0) :+ 1
    ib = ipwadh_cnt_le(T, t1) :+ 1
    xbe = J(n, p, 0)
    xbe[evi, .] = xbar[jj, .]
    s = (w :* ev) :* (X - xbe) - r :* (X :* (Lam[ib] - Lam[ia]) - (Gam[ib, .] - Gam[ia, .]))
}

void ipwadh_se_start()
{
    external transmorphic IPWSE
    IPWSE = asarray_create()
    asarray(IPWSE, "nfit", 0)
}

void ipwadh_se_clear()
{
    external transmorphic IPWSE
    IPWSE = asarray_create()
}

// row-level data (the reduced data set: one observation per usable original row). vparts = "E0 J0 dA Ec0 dC" (the stored parts)
void ipwadh_se_rows(string scalar vid, string scalar vt0, string scalar vlen, string scalar vadh, string scalar vorig,
                    string scalar vchg, string scalar vparts, real scalar ncens)
{
    external transmorphic IPWSE
    real colvector id, t0, len, adh, orig, chg, p
    real matrix P

    id   = st_data(., vid)
    t0   = st_data(., vt0)
    len  = st_data(., vlen)
    adh  = st_data(., vadh)
    orig = st_data(., vorig)
    chg  = st_data(., vchg)
    P    = st_data(., tokens(vparts))
    p = order((id, t0), (1, 2))
    asarray(IPWSE, "p", p)
    id = id[p]
    asarray(IPWSE, "id", id)
    asarray(IPWSE, "t0", t0[p])
    asarray(IPWSE, "len", len[p])
    asarray(IPWSE, "adh", adh[p])
    asarray(IPWSE, "orig", orig[p])
    asarray(IPWSE, "chg", chg[p])
    asarray(IPWSE, "P", P[p, .])
    asarray(IPWSE, "info", panelsetup(id, 1))
    asarray(IPWSE, "ncens", ncens)
}

// one Poisson rate model: fam 1 = adherence change, 2 = censoring cause; role 1 = numerator (reduced), 2 = denominator (full);
// cause = index of the censoring cause (1..ncens) for fam 2, 0 for fam 1
void ipwadh_se_addfit(real scalar f, real scalar fam, real scalar role, real scalar cause, string scalar zlist,
                      string scalar vmask, string scalar vy, string scalar bm, string scalar omm, string scalar Vm)
{
    external transmorphic IPWSE
    real colvector p, mask, y, keep, th
    real rowvector b, om
    real matrix Z, V, Vfull
    real scalar k, kk

    p = asarray(IPWSE, "p")
    mask = st_data(., vmask)
    y = st_data(., vy)
    if (zlist != "") Z = st_data(., tokens(zlist))
    else Z = J(rows(mask), 0, 0)
    Z = (Z, J(rows(mask), 1, 1))
    b = st_matrix(bm)
    om = st_matrix(omm)
    Vfull = st_matrix(Vm)
    keep = selectindex(om' :== 0)
    th = (b[1, keep'])'
    V = Vfull[keep, keep]
    if (rows(th) != cols(Z)) {
        errprintf("stcox_ipwadh internal: fit %f has %f coefficients but %f design columns\n", f, rows(th), cols(Z))
        exit(498)
    }
    asarray(IPWSE, "Z" + strofreal(f), Z[p, .])
    asarray(IPWSE, "mask" + strofreal(f), mask[p])
    asarray(IPWSE, "y" + strofreal(f), y[p])
    asarray(IPWSE, "th" + strofreal(f), th)
    asarray(IPWSE, "V" + strofreal(f), V)
    asarray(IPWSE, "fam" + strofreal(f), fam)
    asarray(IPWSE, "role" + strofreal(f), role)
    asarray(IPWSE, "cause" + strofreal(f), cause)
    if (f > asarray(IPWSE, "nfit")) asarray(IPWSE, "nfit", f)
}

// rates of all rows from the refitted models, and the five weight parts recomputed from them (E0 J0 dA Ec0 dC)
real matrix ipwadh_se_parts()
{
    external transmorphic IPWSE
    real scalar f, nfit, n, c, nc, i, a, b, fam, role, cs
    real colvector id, len, chg, lan, lad, lam, mask, e1, j1, ce, cj, cjl, edt, dct, ced
    real matrix info, MN, MD, Z

    nfit = asarray(IPWSE, "nfit")
    nc = asarray(IPWSE, "ncens")
    len = asarray(IPWSE, "len")
    chg = asarray(IPWSE, "chg")
    info = asarray(IPWSE, "info")
    n = rows(len)
    lan = J(n, 1, 0)
    lad = J(n, 1, 0)
    MN = J(n, max((nc, 1)), 0)
    MD = J(n, max((nc, 1)), 0)
    for (f = 1; f <= nfit; f++) {
        Z = asarray(IPWSE, "Z" + strofreal(f))
        mask = asarray(IPWSE, "mask" + strofreal(f))
        fam = asarray(IPWSE, "fam" + strofreal(f))
        role = asarray(IPWSE, "role" + strofreal(f))
        cs = asarray(IPWSE, "cause" + strofreal(f))
        lam = mask :* exp(Z * asarray(IPWSE, "th" + strofreal(f)))
        if (fam == 1) {
            if (role == 1) lan = lan + lam
            else lad = lad + lam
        }
        else {
            if (role == 1) MN[., cs] = MN[., cs] + lam
            else MD[., cs] = MD[., cs] + lam
        }
    }
    if (min(lan) <= 0 | min(lad) <= 0) {
        errprintf("stcox_ipwadh: a row is not covered by an adherence model (non-positive rate)\n")
        exit(498)
    }
    asarray(IPWSE, "lan", lan)
    asarray(IPWSE, "lad", lad)
    asarray(IPWSE, "MN", MN)
    asarray(IPWSE, "MD", MD)
    e1 = -(lan - lad) :* len
    j1 = (chg :== 1) :* ln(lan :/ lad)
    edt = J(n, 1, 0)
    dct = J(n, 1, 0)
    for (c = 1; c <= nc; c++) {
        edt = edt - (MN[., c] - MD[., c]) :* len
        dct = dct + (MN[., c] - MD[., c])
    }
    ce = J(n, 1, 0)
    cj = J(n, 1, 0)
    cjl = J(n, 1, 0)
    ced = J(n, 1, 0)
    for (i = 1; i <= rows(info); i++) {
        a = info[i, 1]
        b = info[i, 2]
        ce[|a \ b|] = runningsum(e1[|a \ b|])
        cj[|a \ b|] = runningsum(j1[|a \ b|])
        ced[|a \ b|] = runningsum(edt[|a \ b|])
        if (b > a) cjl[|a + 1 \ b|] = cj[|a \ b - 1|]
    }
    return((ce - e1, cjl, lan - lad, ced - edt, dct))
}

// the main computation. All piece-level variables are the full data set (restored); touse = estimation sample of stcox.
//   evalat 1 = end, 2 = mid, 3 = start ; beta / omit / Vcox : e(b), its omit row, e(V) of the stcox run (cluster-robust)
void ipwadh_se_final(string scalar vt0, string scalar vt1, string scalar vd, string scalar vw, string scalar vpsw,
                     string scalar vorig, string scalar vstr, string scalar vtouse, string scalar vtt0, string scalar vtt1,
                     string scalar vx, string scalar bm, string scalar omm, string scalar Vcox, string scalar Vout,
                     real scalar evalat, real scalar fpc)
{
    external transmorphic IPWSE
    real scalar n, nr, pb, nid, i, a, b, f, nfit, fam, role, cs, sgn, G, fac, maxdP, maxdw, maxdV, maxU, sidx, Gs
    real colvector t0, t1, d, w, psw, orig, cstr, tt0, tt1, beta, keep, id, len, chg, rt0, ko, ks, pos, kk, tev, x, g, lw, wchk
    real colvector o, pidx, pp, pio, kko, strs, ix, rowid, ordr, mask, y, mu, resid, lam, lan, lad
    real matrix X, s, A, Aall, Ainv, Ui, psi, corrm, info, infop, Rk, Cur, Fut, Q, P, Pn, kap, Bf, Z, V, Sf, MN, MD, gs, Ve, Vsel, Vnaive, Vadj, Vfull, sb, Ab
    real colvector U, Uall, th
    real rowvector om, bb
    real matrix xs

    // ---------- data ----------
    t0 = st_data(., vt0, vtouse)
    t1 = st_data(., vt1, vtouse)
    d = st_data(., vd, vtouse)
    w = st_data(., vw, vtouse)
    psw = st_data(., vpsw, vtouse)
    orig = st_data(., vorig, vtouse)
    cstr = st_data(., vstr, vtouse)
    tt0 = st_data(., vtt0, vtouse)
    tt1 = st_data(., vtt1, vtouse)
    X = st_data(., tokens(vx), vtouse)
    n = rows(t0)
    pb = cols(X)
    om = st_matrix(omm)
    keep = selectindex(om' :== 0)
    bb = st_matrix(bm)
    beta = (bb[1, keep'])'
    if (rows(beta) != pb) {
        errprintf("stcox_ipwadh internal: %f retained coefficients but %f design columns\n", rows(beta), pb)
        exit(498)
    }
    if (max(abs(t0 - tt0)) > 1e-9 * max((1, max(abs(t1)))) | max(abs(t1 - tt1)) > 1e-9 * max((1, max(abs(t1))))) {
        errprintf("stcox_ipwadh: the analysis times _t0/_t differ from the time0()/time1() variables (stset origin() or scale()?); not supported\n")
        exit(498)
    }

    // ---------- rows: refitted models, recomputed parts, comparison with the stored parts ----------
    P = ipwadh_se_parts()
    Pn = asarray(IPWSE, "P")
    maxdP = max(abs(P - Pn) :/ (1 :+ abs(Pn)))
    id = asarray(IPWSE, "id")
    len = asarray(IPWSE, "len")
    chg = asarray(IPWSE, "chg")
    rt0 = asarray(IPWSE, "t0")
    info = asarray(IPWSE, "info")
    nr = rows(id)
    nid = rows(info)
    MN = asarray(IPWSE, "MN")
    MD = asarray(IPWSE, "MD")
    lan = asarray(IPWSE, "lan")
    lad = asarray(IPWSE, "lad")

    // ---------- map the pieces of the Cox sample to their rows ----------
    ko = order(asarray(IPWSE, "orig"), 1)
    ks = asarray(IPWSE, "orig")
    ks = ks[ko]
    pos = ipwadh_cnt_lt(ks, orig) :+ 1
    if (max(pos) > nr) {
        errprintf("stcox_ipwadh: pieces of the Cox sample without a usable original row (%s)\n", vorig)
        exit(498)
    }
    if (any(ks[pos] :!= orig)) {
        errprintf("stcox_ipwadh: pieces of the Cox sample without a usable original row (%s)\n", vorig)
        exit(498)
    }
    kk = ko[pos]
    tev = (evalat == 1 ? tt1 : (evalat == 2 ? (tt0 + tt1) :/ 2 : tt0))
    x = tev - rt0[kk]
    lw = P[kk, 1] - P[kk, 3] :* x + P[kk, 2] + P[kk, 4] - P[kk, 5] :* x
    wchk = exp(lw)
    maxdw = max(abs(wchk - psw) :/ psw)
    g = (abs(w - psw) :<= 1e-9 :* psw)

    // ---------- Cox quantities (Breslow, pweights), by Cox stratum ----------
    s = J(n, pb, 0)
    Aall = J(pb, pb, 0)
    Uall = J(pb, 1, 0)
    strs = uniqrows(cstr)
    for (sidx = 1; sidx <= rows(strs); sidx++) {
        ix = selectindex(cstr :== strs[sidx])
        ipwadh_cox_stratum(t0[ix], t1[ix], d[ix], w[ix], X[ix, .], beta, sb, Ab, U)
        s[ix, .] = sb
        Aall = Aall + Ab
        Uall = Uall + U
    }
    maxU = max(abs(Uall))
    Ainv = invsym(Aall)

    // ---------- person-level Cox scores U_i ----------
    pidx = J(nr, 1, 0)
    for (i = 1; i <= nid; i++) pidx[|info[i, 1] \ info[i, 2]|] = J(info[i, 2] - info[i, 1] + 1, 1, i)
    pp = pidx[kk]
    o = order(pp, 1)
    pio = pp[o]
    infop = panelsetup(pio, 1)
    Gs = rows(infop)
    Ui = J(nid, pb, 0)
    pio = pp[o]
    Ui[pio[infop[., 1]], .] = panelsum(s[o, .], infop)
    G = Gs
    fac = (fpc ? G / (G - 1) : 1)
    Vnaive = fac :* Ainv * (Ui' * Ui) * Ainv

    // ---------- row-level sums of the score contributions: current row (weighted by x), future rows ----------
    gs = s :* g
    Rk = J(nr, pb, 0)
    Cur = J(nr, pb, 0)
    o = order(kk, 1)
    infop = panelsetup(kk[o], 1)
    kko = kk[o]
    rowid = kko[infop[., 1]]
    Rk[rowid, .] = panelsum(gs[o, .], infop)
    Cur[rowid, .] = panelsum(gs[o, .] :* x[o], infop)
    Fut = J(nr, pb, 0)
    for (i = 1; i <= nid; i++) {
        a = info[i, 1]
        b = info[i, 2]
        Fut[|a, 1 \ b, pb|] = J(b - a + 1, 1, 1) * colsum(Rk[|a, 1 \ b, pb|]) - ipwadh_cs(Rk[|a, 1 \ b, pb|])
    }
    Q = len :* Fut + Cur

    // ---------- B_f V_f S_if over the rate models ----------
    nfit = asarray(IPWSE, "nfit")
    corrm = J(nid, pb, 0)
    for (f = 1; f <= nfit; f++) {
        Z = asarray(IPWSE, "Z" + strofreal(f))
        mask = asarray(IPWSE, "mask" + strofreal(f))
        y = asarray(IPWSE, "y" + strofreal(f))
        th = asarray(IPWSE, "th" + strofreal(f))
        V = asarray(IPWSE, "V" + strofreal(f))
        fam = asarray(IPWSE, "fam" + strofreal(f))
        role = asarray(IPWSE, "role" + strofreal(f))
        cs = asarray(IPWSE, "cause" + strofreal(f))
        sgn = (role == 1 ? -1 : 1)
        if (fam == 1) {
            lam = (role == 1 ? lan : lad)
            kap = sgn :* lam :* Q - sgn :* chg :* Fut
        }
        else {
            lam = (role == 1 ? MN[., cs] : MD[., cs])
            kap = sgn :* lam :* Q
        }
        kap = mask :* kap
        Bf = kap' * Z
        mu = mask :* exp(Z * th) :* len
        resid = mask :* (y - mu)
        Sf = panelsum(Z :* resid, info)
        corrm = corrm + (Sf * V) * Bf'
    }
    psi = Ui + corrm
    Vadj = fac :* Ainv * (psi' * psi) * Ainv

    // ---------- comparison with the cluster-robust e(V) of stcox ----------
    Ve = st_matrix(Vcox)
    Vsel = Ve[keep, keep]
    maxdV = max(abs(diagonal(Vnaive) - diagonal(Vsel)) :/ diagonal(Vsel))

    Vfull = J(rows(Ve), cols(Ve), 0)
    Vfull[keep, keep] = Vadj
    st_matrix(Vout, Vfull)
    st_numscalar("__ipwse_maxdP", maxdP)
    st_numscalar("__ipwse_maxdw", maxdw)
    st_numscalar("__ipwse_maxdV", maxdV)
    st_numscalar("__ipwse_maxU", maxU)
    st_numscalar("__ipwse_G", G)
    st_numscalar("__ipwse_nclip", sum(g :== 0))
    st_numscalar("__ipwse_npieces", n)
}

end
