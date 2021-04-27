*! version 1.0.1  27apr2021  Ben Jann
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
    transmorphic   ls
}

struct mm_areg_struct scalar mm_areg(
    real colvector   y,
    real colvector   id,
    | real matrix    X,
      real colvector w,
      real scalar    sort,
      real scalar    qd) 
{
    struct mm_areg_struct scalar t
    real colvector ym, yd
    real matrix    Xm, Xd
    pragma unset   ym
    pragma unset   Xm
    
    // setup
    if (args()<3) X = J(rows(y), 0, .)
    if (args()<4) w = 1
    t.N = t.rss = t.s = t.r2 = t.XXinv = t.V = t.ymean = t.means = .z
    t.y = &y
    if (X==.) t.X = &(J(rows(y), 0, .))
    else      t.X = &X
    t.w = &w
    t.k = cols(*t.X)
    
    // estimate
    t.levels = _mm_areg_m(y, *t.X, w, id, sort, ym, Xm, t.n)
    yd = y - ym; Xd = *t.X - Xm
    t.ls = mm_ls(yd, Xd, w, 0, qd, 0)
    yd = J(0, 1, .); Xd = J(0, 0, .) // no longer needed; free memory
    
    // recover constant and fixed-effects
    t.a = mm_areg_ymean(t)
    if (t.k) t.a = t.a - mm_ls_xb(t.ls, mm_areg_means(t))
    t.u = ym :- t.a
    if (t.k) t.u = t.u - mm_ls_xb(t.ls, Xm)
    
    // return
    return(t)
}

real colvector _mm_areg_m(real colvector y, real matrix X, real colvector w,
    real colvector id, real scalar sort, real colvector ym, real matrix Xm,
    real colvector n)
{
    real matrix    info
    real scalar    j
    real colvector p, levels
    real colvector ws
    
    j = cols(X)
    Xm = J(rows(X), j, .)
    if (sort) {
        p = order(id, 1)
        n = _mm_panels(id[p])
        info = mm_colrunsum(n)
        levels = (id[p])[info] 
        ws = (rows(w)==1 ? w : w[p])
        ym = J(rows(y),1,.)
        ym[p] = __mm_areg_m(y[p], ws, info)
        for (; j; j--) Xm[p,j] = __mm_areg_m(X[p,j], ws, info)
        return(levels)
    }
    n = _mm_panels(id)
    info = mm_colrunsum(n)
    levels = id[info]
    ym = __mm_areg_m(y, w, info)
    for (; j; j--) Xm[,j] = __mm_areg_m(X[,j], w, info)
    return(levels)
}

real colvector __mm_areg_m(real colvector y, real colvector w, real matrix info)
{
    real scalar    i, n, a, b, ww
    real colvector m

    m = y
    if (rows(m)<1) return(m)
    ww = (rows(w)!=1)
    n = rows(info)
    b = 0
    for (i=1; i<=n; i++) {
        a = b + 1
        b = info[i]
        m[|a \ b|] = J(b-a+1, 1, mean(m[|a \ b|], ww ? w[|a \ b|] : w))
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
