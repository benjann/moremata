*! version 1.0.4  15sep2021  Ben Jann
version 9.2
mata:

// Main functions

real colvector mm_aregfit(real colvector y, real colvector id, | real matrix X,
    real colvector w, real scalar sort, real scalar qd)
{
    if (args()<3) X = .
    if (args()<4) w = 1
    return(mm_areg_b(mm_areg(y, id, X, w, sort, qd)))
}

struct mm_areg_struct {
    pointer(real colvector) scalar y, w
    pointer(real matrix) scalar    X
    real scalar                    qd, k, a, N, ymean, rss, s, r2
    real colvector                 yd, xb, e, u
    real rowvector                 means
    real matrix                    Xd, XXinv, V
    struct mm_ls_struct scalar     ls
    pointer(struct mm_areg_struct_g scalar) scalar g
}

struct mm_areg_struct_g {
    real scalar    sort             // 1 if data requires sorting, 0 else
    real colvector p                // sort index; only filled in if sort==1
    real colvector idx, levels, n
}

struct mm_areg_struct scalar mm_areg(
    real colvector   y,
    real colvector   id,
    | real matrix    X,
      real colvector w,
      real scalar    sort,
      real scalar    qd)
{
    if (args()<3) X = .
    if (args()<4) w = 1
    if (rows(id)!=rows(y)) _error(3200)
    return(_mm_areg(_mm_areg_g(id, sort), y, X, w, qd))
}

struct mm_areg_struct_g scalar _mm_areg_g(real colvector id, real scalar sort)
{
    struct mm_areg_struct_g scalar g
    
    g.sort = (sort!=0)
    if (g.sort) {
        g.p = mm_order(id, 1, 1) // stable sort order
        __mm_areg_g(g, id[g.p])
    }
    else {
        __mm_areg_g(g, id)
    }
    return(g)
}

void __mm_areg_g(struct mm_areg_struct_g scalar g, real colvector id)
{
    real scalar r
    
    r = rows(id)
    if (r==0) {
        g.idx = g.levels = g.n = J(0,1,.)
    }
    else if (r==1) {
        g.idx = 1
        g.levels = id
        g.n = r
    }
    else {
        g.idx = (id :!= (id[|2\.|] \ id[1]))
        g.idx[r] = 1 // should id be constant
        g.idx = select(1::r, g.idx) // index of last obs in each group
        g.levels = id[g.idx]
        g.n = g.idx - (0 \ g.idx)[|1 \ rows(g.idx)|] // group sizes
    }
}

struct mm_areg_struct scalar _mm_areg(
    struct mm_areg_struct_g scalar g,
    real colvector                 y,
    | real matrix                  X,
      real colvector               w,
      real scalar                  qd)
{
    struct mm_areg_struct scalar t
    
    // setup
    if (args()<3) X = .
    if (args()<4) w = 1
    t.qd = (qd!=0)
    t.a = t.xb = t.e = t.u = .z
    t.N = t.rss = t.s = t.r2 = t.XXinv = t.V = t.ymean = t.means = .z
    t.g = &g
    t.y = &y
    if (X==.) t.X = &(J(rows(y), 0, .))
    else      t.X = &X
    t.w = &w
    t.k = cols(*t.X)
    
    // transform data
    _mm_areg_demean(t.yd, t.Xd, y, *t.X, w, g, t.qd)
    
    // estimate
    t.ls = mm_ls(t.yd, t.Xd, w, 0, t.qd)
    
    // return
    return(t)
}

void _mm_areg_demean(real colvector yd, real matrix Xd, 
    real colvector y, real matrix X, real colvector w,
    struct mm_areg_struct_g scalar g, | real scalar qd)
{
    real matrix M
    
    M = (y,X)
    _mm_areg_gmean(g, M, w, qd)
    yd = y - M[,1]
    if (cols(X)) Xd = X - M[,2..cols(M)]
    else         Xd = X
}

void _mm_areg_gmean(struct mm_areg_struct_g scalar g, real matrix X,
    real colvector w, | real scalar qd)
{
    if (g.sort) {
        _collate(X, g.p)
        if (qd) __mm_areg_gmean(X, rows(w)==1 ? w : w[g.p], g.idx)
        else    __mm_areg_gmean_dbl(X, rows(w)==1 ? w : w[g.p], g.idx)
        _collate(X, invorder(g.p))
    }
    else {
        if (qd) __mm_areg_gmean(X, w, g.idx)
        else    __mm_areg_gmean_dbl(X, w, g.idx)
    }
}

void __mm_areg_gmean(real matrix X, real colvector w, real colvector idx)
{
    real scalar    i, n, a, b, k, W
    real colvector ww
    real matrix    ab
    
    if (rows(X)<1) return
    n = rows(idx)
    b = 0
    if (rows(w)==1) {
        for (i=1; i<=n; i++) {
            a = b + 1
            b = idx[i]
            k = b - a
            if (!k) continue // no averaging necessary if less than two obs
            k++
            ab = (a,1 \ b,.)
            X[|ab|] = J(k, 1, quadcolsum(X[|ab|])/k)
        }
        return
    }
    for (i=1; i<=n; i++) {
        a = b + 1
        b = idx[i]
        k = b - a
        if (!k) continue // no averaging necessary if less than two obs
        k++
        ab = (a,1 \ b,.)
        ww = w[|ab|]
        W  = quadsum(ww)
        if (W) X[|ab|] = J(k, 1, quadcross(ww/W, X[|ab|]))
        else   X[|ab|] = J(k, 1, quadcolsum(X[|ab|])/k)
            // using unweighted average if sum of weights is 0
    }
}

void __mm_areg_gmean_dbl(real matrix X, real colvector w, real colvector idx)
{   // double-precision variant of __mm_areg_gmean()
    real scalar    i, n, a, b, k, W
    real colvector ww
    real matrix    ab

    if (rows(X)<1) return
    n = rows(idx)
    b = 0
    if (rows(w)==1) {
        for (i=1; i<=n; i++) {
            a = b + 1
            b = idx[i]
            k = b - a
            if (!k) continue // no averaging necessary if less than two obs
            k++
            ab = (a,1 \ b,.)
            X[|ab|] = J(k, 1, colsum(X[|ab|])/k)
        }
        return
    }
    for (i=1; i<=n; i++) {
        a = b + 1
        b = idx[i]
        k = b - a
        if (!k) continue // no averaging necessary if less than two obs
        k++
        ab = (a,1 \ b,.)
        ww = w[|ab|]
        W  = sum(ww)
        if (W) X[|ab|] = J(k, 1, cross(ww/W, X[|ab|]))
        else   X[|ab|] = J(k, 1, colsum(X[|ab|])/k)
            // using unweighted average if sum of weights is 0
    }
}

// Results

real colvector mm_areg_b(struct mm_areg_struct scalar t)
{
    return(mm_areg_beta(t) \ mm_areg_alpha(t))
}

real colvector mm_areg_beta(struct mm_areg_struct scalar t)
{
    return(mm_ls_b(t.ls))
}

real scalar mm_areg_alpha(struct mm_areg_struct scalar t)
{
    if (t.a==.z) {
        if (t.k) t.a = mm_areg_ymean(t) - mm_ls_xb(t.ls, mm_areg_means(t))
        else     t.a = mm_areg_ymean(t)
    }
    return(t.a)
}

real colvector mm_areg_xb(struct mm_areg_struct scalar t, | real matrix X)
{
    if (args()==2) {
        if (cols(X)!=t.k) _error(3200)
        if (t.k) return(mm_ls_xb(t.ls, X) :+ mm_areg_alpha(t))
        return(J(rows(X), 1, mm_areg_alpha(t)))
    }
    if (t.xb==.z) {
        if (t.k) t.xb = mm_ls_xb(t.ls, *t.X) :+ mm_areg_alpha(t)
        else     t.xb = J(rows(*t.y), 1, mm_areg_alpha(t))
    }
    return(t.xb)
}

real colvector mm_areg_ue(struct mm_areg_struct scalar t)
{
    return(*t.y - mm_areg_xb(t))
}

real colvector mm_areg_xbu(struct mm_areg_struct scalar t)
{
    return(mm_areg_xb(t) + mm_areg_u(t))
}

real colvector mm_areg_u(struct mm_areg_struct scalar t)
{
    if (t.u==.z) {
        t.u = *t.y - mm_areg_xb(t)
        _mm_areg_gmean(*t.g, t.u, *t.w, t.qd)
        // could also computed t.u = *t.y - mm_areg_e(t) - mm_areg_xb(t), but
        // would still have to take care of roundoff error to make sure that
        // t.u is truly constant within group
    }
    return(t.u)
}

real colvector mm_areg_e(struct mm_areg_struct scalar t)
{
    if (t.e==.z) {
        if (t.k) t.e = t.yd - mm_ls_xb(t.ls)
        else     t.e = t.yd
    }
    return(t.e)
}

real scalar mm_areg_s(struct mm_areg_struct scalar t)
{
    if (t.s==.z) {
        t.s = sqrt(mm_areg_rss(t) * 
            editmissing(1 / (mm_areg_N(t) - (t.k - mm_ls_k_omit(t.ls))
                 - rows(t.g->idx)), 0))
    }
    return(t.s)
}

real scalar mm_areg_r2(struct mm_areg_struct scalar t)
{
    if (t.r2==.z) {
        t.r2 = 1 - mm_areg_rss(t) / 
                   quadcross(*t.w, (*t.y :- mm_areg_ymean(t)):^2)
    }
    return(t.r2)
}

real colvector mm_areg_se(struct mm_areg_struct scalar t)
{
    return(sqrt(diagonal(mm_areg_V(t))))
}

real matrix mm_areg_V(struct mm_areg_struct scalar t)
{
    if (t.V==.z) t.V = mm_areg_XXinv(t)*mm_areg_s(t)^2
    return(t.V)
}

real matrix mm_areg_XXinv(struct mm_areg_struct scalar t)
{
    if (t.XXinv==.z) {
        t.XXinv = mm_ls_XXinv(t.ls) \ -(mm_areg_means(t) * mm_ls_XXinv(t.ls))
        t.XXinv = (t.XXinv, (t.XXinv[t.k+1,]' \
                   1/mm_areg_N(t) - mm_areg_means(t) * t.XXinv[t.k+1,]'))
    }
    return(t.XXinv)
}

real scalar mm_areg_rss(struct mm_areg_struct scalar t)
{
    if (t.rss==.z) {
        t.rss = quadcross(*t.w, mm_areg_e(t):^2)
    }
    return(t.rss)
}

real scalar mm_areg_ymean(struct mm_areg_struct scalar t)
{
    real scalar n
    
    if (t.ymean==.z) {
        n = rows(*t.y)
        if (n==0)               t.ymean = .
        else if (rows(*t.w)==1) t.ymean = quadsum(*t.y) / n
        else                    t.ymean = quadsum(*t.w :* *t.y) / mm_areg_N(t)
    }
    return(t.ymean)
}

real rowvector mm_areg_means(struct mm_areg_struct scalar t)
{
    if (t.means==.z) {
        if      (t.k==0)        t.means = J(1, 0, .)
        else if (rows(*t.X)==0) t.means = J(1, t.k, .)
        else                    t.means = quadcross(*t.w, *t.X) / mm_areg_N(t)
    }
    return(t.means)
}

real matrix mm_areg_yd(struct mm_areg_struct scalar t)
{
    return(t.yd)
}

real matrix mm_areg_Xd(struct mm_areg_struct scalar t)
{
    return(t.Xd)
}

real colvector mm_areg_omit(struct mm_areg_struct scalar t)
{
    return(mm_ls_omit(t.ls) \ 0)
}

real colvector mm_areg_k_omit(struct mm_areg_struct scalar t)
{
    return(mm_ls_k_omit(t.ls))
}

real scalar mm_areg_N(struct mm_areg_struct scalar t)
{
    if (t.N==.z) {
        t.N = rows(*t.w)==1 ? *t.w * rows(*t.y) : quadsum(*t.w)
    }
    return(t.N)
}

real colvector mm_areg_levels(struct mm_areg_struct scalar t)
{
    return(t.g->levels)
}

real colvector mm_areg_k_levels(struct mm_areg_struct scalar t)
{
    return(rows(t.g->idx))
}

real colvector mm_areg_n(struct mm_areg_struct scalar t)
{
    return(t.g->n)
}

end
