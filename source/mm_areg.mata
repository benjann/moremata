*! version 1.0.0  23apr2021  Ben Jann
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
    pointer(real colvector) scalar y
    pointer(real matrix) scalar    X
    real scalar                    a, N, ymean, rss, s, r2
    real colvector                 u, levels, n
    real rowvector                 means
    real matrix                    XXinv, V
    struct mm_ls_struct scalar     ls
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
    real colvector ym
    real matrix    Xm
    pragma unset   ym
    pragma unset   Xm
    
    // setup
    if (args()<3) X = J(rows(y), 0, .)
    if (args()<4) w = 1
    t.N = t.rss = t.s = t.r2 = t.XXinv = t.V = t.ymean = t.means = .z
    t.y = &y
    if (X==.) t.X = &(J(rows(y), 0, .))
    else      t.X = &X
    
    // estimate
    t.levels = _mm_areg_m(y, *t.X, w, id, sort, ym, Xm, t.n)
    t.ls = mm_ls(y-ym, *t.X-Xm, w, 0, qd, 0)
    t.ls.y = t.ls.X = NULL // no longer needed; free memory
    
    // recover constant and fixed-effects
    t.a = mm_areg_ymean(t) - (t.ls.k ? mm_ls_xb(t.ls, mm_areg_means(t)) : 0)
    t.u = ym - (t.a :+ (t.ls.k ? mm_ls_xb(t.ls, Xm) : J(rows(*t.y),1,0)))
    
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
        else               t.ymean = quadcross(*t.ls.w, *t.y) / mm_areg_N(t)
    }
    return(t.ymean)
}

real rowvector mm_areg_means(struct mm_areg_struct scalar t)
{
    if (t.means==.z) {
        if      (t.ls.k==0)     t.means = J(1, 0, .)
        else if (rows(*t.X)==0) t.means = J(1, t.ls.k, .)
        else                    t.means = quadcross(*t.ls.w, *t.X) / mm_areg_N(t)
    }
    return(t.means)
}

real colvector mm_areg_b(struct mm_areg_struct scalar t)
{
    return(t.ls.b \ t.a)
}

real colvector mm_areg_xb(struct mm_areg_struct scalar t, | real matrix X)
{
    if (args()==2) {
        return((t.ls.k ? mm_ls_xb(t.ls, X) : J(rows(*t.y),1,0)) :+ t.a)
    }
    return((t.ls.k ? mm_ls_xb(t.ls, *t.X) : J(rows(*t.y),1,0)) :+ t.a)
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
        t.N = rows(*t.ls.w)==1 ? *t.ls.w * rows(*t.y) : quadsum(*t.ls.w)
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
        t.rss = quadcross(*t.ls.w, (*t.y - mm_areg_xbu(t)):^2)
    }
    return(t.rss)
}

real scalar mm_areg_s(struct mm_areg_struct scalar t)
{
    if (t.s==.z) {
        t.s = sqrt(mm_areg_rss(t) * 
            editmissing(1 / (mm_areg_N(t)-t.ls.kadj-rows(t.levels)), 0))
    }
    return(t.s)
}

real scalar mm_areg_r2(struct mm_areg_struct scalar t)
{
    if (t.r2==.z) {
        t.r2 = 1 - mm_areg_rss(t) / 
                   quadcross(*t.ls.w, (*t.y :- mm_areg_ymean(t)):^2)
    }
    return(t.r2)
}

real matrix mm_areg_XXinv(struct mm_areg_struct scalar t)
{
    if (t.XXinv==.z) {
        t.XXinv = t.ls.XXinv \ -(mm_areg_means(t) * t.ls.XXinv)
        t.XXinv = (t.XXinv, (t.XXinv[t.ls.k+1,]' \
                   1/mm_areg_N(t) - mm_areg_means(t) * t.XXinv[t.ls.k+1,]'))
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
