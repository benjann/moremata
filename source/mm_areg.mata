*! version 1.0.3  06sep2021  Ben Jann
version 9.2
mata:

real colvector mm_aregfit(real colvector y, real colvector id, | real matrix X,
    real colvector w, real scalar sort, real scalar qd)
{
    if (args()<3) X = J(rows(y), 0, .)
    if (args()<4) w = 1
    return(mm_areg_b(mm_areg(y, id, X, w, sort, qd)))
}

struct mm_areg_struct {
    pointer(real colvector) scalar y, w
    pointer(real matrix) scalar    X
    real scalar    k, a, N, ymean, rss, s, r2
    real colvector u, levels, n
    real rowvector means
    real matrix    XXinv, V
    struct mm_ls_struct scalar ls
}

struct mm_areg_struct_grps {
    real colvector p
    real colvector idx, levels, n
}

struct mm_areg_struct_data {
    real colvector ym, yd
    real matrix    Xm, Xd
}

struct mm_areg_struct scalar mm_areg(
    real colvector   y,
    real colvector   id,
    | real matrix    X,
      real colvector w,
      real scalar    sort,
      real scalar    qd,
      struct mm_areg_struct_grps g,
      // - if empty g is provided, g will be filled in with info on groups
      // - if filled-in g is provided, info on groups will be taken from g
      // - g is considered empty if rows(g.idx)==0
      struct mm_areg_struct_data d
      // - if empty d is provided, d will be filled in with transformed data
      // - if filled-in d is provided, transformed data will be taken from d
      // - status of d will be evaluated separately for y and X
      // - d will be considered empty for y, if rows(d.ym)==0
      // - d will be considered empty for X, if rows(d.Xm)==0
      )
{
    struct mm_areg_struct scalar t
    
    // setup
    if (args()<3) X = J(rows(y), 0, .)
    if (args()<4) w = 1
    if (length(g)==0) g = mm_areg_struct_grps(1)
    if (length(d)==0) d = mm_areg_struct_data(1)
    t.N = t.rss = t.s = t.r2 = t.XXinv = t.V = t.ymean = t.means = .z
    t.y = &y
    if (X==.) t.X = &(J(rows(y), 0, .))
    else      t.X = &X
    t.w = &w
    t.k = cols(*t.X)
    if (rows(id)!=rows(y)) _error(3200)
    
    // collect information on groups
    if (rows(g.idx)==0) _mm_areg_grps(g, id, sort)
    t.levels = g.levels; t.n = g.n
    
    // transform data
    if (rows(d.ym)==0) {
        d.ym = _mm_areg_gmean(g, y, w, sort)
        d.yd = y - d.ym
    }
    if (rows(d.Xm)==0) {
        d.Xm = _mm_areg_gmean(g, *t.X, w, sort)
        d.Xd = *t.X - d.Xm
    }
    
    // estimate
    t.ls = mm_ls(d.yd, d.Xd, w, 0, qd, 0)
    t.ls.y = t.ls.X = t.ls.w = NULL // release pointers to data (save memory)
    
    // recover constant and fixed-effects
    t.a = mm_areg_ymean(t)
    if (t.k) t.a = t.a - mm_ls_xb(t.ls, mm_areg_means(t))
    t.u = d.ym :- t.a
    if (t.k) t.u = t.u - mm_ls_xb(t.ls, d.Xm)
    
    // return
    return(t)
}

void _mm_areg_grps(struct mm_areg_struct_grps scalar g, real colvector id, 
    real scalar sort)
{
    if (sort) __mm_areg_grps(g, id[g.p = order(id, 1)])
    else      __mm_areg_grps(g, id)
}

void __mm_areg_grps(struct mm_areg_struct_grps scalar g, real colvector id)
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

real matrix _mm_areg_gmean(struct mm_areg_struct_grps scalar g, 
    real matrix X, real colvector w, real scalar sort)
{
    real scalar    j
    real colvector ws
    real matrix    Xm
    
    j = cols(X)
    if (j==0) return(J(rows(X), 0, .))
    if (sort) {
        ws = (rows(w)==1 ? w : w[g.p])
        Xm = X[g.p,]
        for (; j; j--) Xm[g.p,j] = __mm_areg_gmean(Xm[,j], ws, g.idx)
    }
    else {
        Xm = J(rows(X), j, .)
        for (; j; j--) Xm[,j] = __mm_areg_gmean(X[,j], w, g.idx)
    }
    return(Xm)
}

real colvector __mm_areg_gmean(real colvector x, real colvector w,
    real colvector idx)
{
    real scalar    i, n, a, b, k, W
    real colvector m, ww

    m = x
    if (rows(m)<1) return(m)
    n = rows(idx)
    b = 0
    if (rows(w)==1) {
        for (i=1; i<=n; i++) {
            a = b + 1
            b = idx[i]
            k = b - a
            if (k) { // no averaging necessary if less than two obs
                k++
                m[|a \ b|] = J(k, 1, quadsum(m[|a \ b|])/k)
            }
        }
        return(m)
    }
    for (i=1; i<=n; i++) {
        a = b + 1
        b = idx[i]
        k = b - a
        if (k) { // no averaging necessary if less than two obs
            k++
            ww = w[|a \ b|]
            W = quadsum(ww)
            if (W) m[|a \ b|] = J(k, 1, quadsum(ww:*m[|a \ b|])/W)
            else   m[|a \ b|] = J(k, 1, quadsum(m[|a \ b|])/k)
                // using unweighted average if sum of weights is 0
        }
    }
    return(m)
}

real scalar mm_areg_ymean(struct mm_areg_struct scalar t)
{
    if (t.ymean==.z) {
        if (rows(*t.y)==0) t.ymean = .
        else               t.ymean = quadcross(*t.w, *t.y) / mm_areg_N(t)
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

real colvector mm_areg_b(struct mm_areg_struct scalar t)
{
    return(mm_ls_b(t.ls) \ t.a)
}

real colvector mm_areg_xb(struct mm_areg_struct scalar t, | real matrix X)
{
    if (args()==2) {
        if (t.k) return(mm_ls_xb(t.ls, X) :+ t.a)
        return(J(rows(*t.y), 1, t.a))
    }
    if (t.k) return(mm_ls_xb(t.ls, *t.X) :+ t.a)
    return(J(rows(*t.y), 1, t.a))
}

real colvector mm_areg_ue(struct mm_areg_struct scalar t)
{
    return(*t.y - mm_areg_xb(t))
}

real colvector mm_areg_xbu(struct mm_areg_struct scalar t)
{
    return(mm_areg_xb(t) + t.u)
}

real colvector mm_areg_u(struct mm_areg_struct scalar t)
{
    return(t.u)
}

real colvector mm_areg_e(struct mm_areg_struct scalar t)
{
    return(*t.y - mm_areg_xbu(t))
}

real colvector mm_areg_levels(struct mm_areg_struct scalar t)
{
    return(t.levels)
}

real colvector mm_areg_k_levels(struct mm_areg_struct scalar t)
{
    return(rows(t.levels))
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

real colvector mm_areg_n(struct mm_areg_struct scalar t)
{
    return(t.n)
}

real scalar mm_areg_rss(struct mm_areg_struct scalar t)
{
    if (t.rss==.z) {
        t.rss = quadcross(*t.w, (*t.y - mm_areg_xbu(t)):^2)
    }
    return(t.rss)
}

real scalar mm_areg_s(struct mm_areg_struct scalar t)
{
    if (t.s==.z) {
        t.s = sqrt(mm_areg_rss(t) * 
            editmissing(1 / (mm_areg_N(t) - (t.k - mm_ls_k_omit(t.ls))
                 - rows(t.levels)), 0))
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

real matrix mm_areg_XXinv(struct mm_areg_struct scalar t)
{
    if (t.XXinv==.z) {
        t.XXinv = mm_ls_XXinv(t.ls) \ -(mm_areg_means(t) * mm_ls_XXinv(t.ls))
        t.XXinv = (t.XXinv, (t.XXinv[t.k+1,]' \
                   1/mm_areg_N(t) - mm_areg_means(t) * t.XXinv[t.k+1,]'))
    }
    return(t.XXinv)
}

real matrix mm_areg_V(struct mm_areg_struct scalar t)
{
    if (t.V==.z) t.V = mm_areg_XXinv(t)*mm_areg_s(t)^2
    return(t.V)
}

real colvector mm_areg_se(struct mm_areg_struct scalar t)
{
    return(sqrt(diagonal(mm_areg_V(t))))
}

end
