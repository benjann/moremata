*! version 1.0.0  27apr2021  Ben Jann
version 9.2
mata:

real matrix mm_aqregfit(real colvector y, real colvector id, | real matrix X,
    real colvector w, real vector tau, real scalar sort, real scalar qd)
{
    if (args()<3) X = J(rows(y), 0, .)
    if (args()<4) w = 1
    if (args()<5) tau = 0.5
    return(mm_aqreg_b(mm_aqreg(y, id, X, w, tau, sort, qd)))
}

struct mm_aqreg_struct {
    pointer(real colvector) scalar y, w
    pointer(real matrix) scalar    X
    real scalar     k, N, ymean
    real colvector  levels, n
    real rowvector  means
    real rowvector  tau, q
    real scalar     beta0, gamma0
    real colvector  alpha, delta
    struct mm_ls_struct scalar L, S
}

struct mm_aqreg_struct_data {
    real colvector ym, yd
    real matrix    Xm, Xd
    real colvector r, rm, rd, e, u
}

struct mm_aqreg_struct scalar mm_aqreg(
    real colvector   y,
    real colvector   id,
    | real matrix    X,
      real colvector w,
      real vector    tau,
      real scalar    sort,
      real scalar    qd,
      struct mm_areg_struct_grps g,
      // - if empty g is provided, g will be filled in with info on groups
      // - if filled-in g is provided, info on groups will be taken from g
      // - g is considered empty if rows(g.idx)==0
      struct mm_aqreg_struct_data d
      // - if empty d is provided, d will be filled in with transformed data
      // - if filled-in d is provided, transformed data will be taken from d
      // - status of d will be evaluated separately for y, X, and r
      // - d will be considered empty for y, if rows(d.ym)==0
      // - d will be considered empty for X, if rows(d.Xm)==0
      // - d will always be considered empty for r, e and u
      )
{
    struct mm_aqreg_struct scalar t
    
    // data
    if (args()<3) X = J(rows(y), 0, .)
    if (args()<4) w = 1
    if (length(g)==0) g = mm_areg_struct_grps(1)
    if (length(d)==0) d = mm_aqreg_struct_data(1)
    t.N = t.ymean = t.means = .z
    t.y = &y
    if (X==.) t.X = &(J(rows(y), 0, .))
    else      t.X = &X
    t.w = &w
    t.k = cols(*t.X)
    if (rows(id)!=rows(y)) _error(3200)
    
    // quantiles
    if (args()<5)    tau = 0.5
    if (rows(tau)>1) t.tau = tau'
    else             t.tau = tau
    if (length(tau)==0) _error(3200)
    if (any(tau:<=0 :| tau:>=1)) _error(3300, "tau must be between 0 and 1")
    
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
    
    // location fit (beta0: global constant; alpha: fixed effects)
    t.L = mm_ls(d.yd, d.Xd, w, 0, qd, 0)
    t.L.y = t.L.X = t.L.w = NULL // release pointers to data (save memory)
    t.beta0 = mm_aqreg_ymean(t)
    if (t.k) t.beta0 = t.beta0 - mm_ls_xb(t.L, mm_aqreg_means(t))
    if (t.k) t.alpha = d.ym - mm_ls_xb(t.L, d.Xm)
    else     t.alpha = d.ym
    t.alpha = t.alpha :- t.beta0
    if (t.k) d.e = d.yd - mm_ls_xb(t.L, d.Xd)
    else     d.e = d.yd
    
    // scale fit (gamma0: global constant; delta: fixed effects)
    d.r = 2 * d.e :* ((d.e:>=0) :- mean(d.e:>=0, w))
    d.rm = _mm_areg_gmean(g, d.r, w, sort)
    d.rd = d.r - d.rm
    t.S = mm_ls(d.rd, d.Xd, w, 0, qd, 0)
    t.S.y = t.S.X = t.S.w = NULL // release pointers to data (save memory)
    t.gamma0 = mean(d.r, w)
    if (t.k) t.gamma0 = t.gamma0 - mm_ls_xb(t.S, mm_aqreg_means(t))
    if (t.k) t.delta = d.rm - mm_ls_xb(t.S, d.Xm)
    else     t.delta = d.rm
    if (t.k) d.u = d.e :/ (mm_ls_xb(t.S, *t.X) + t.delta)
    else     d.u = d.e :/ t.delta
    t.delta = t.delta :- t.gamma0
    
    // quantiles
    if (missing(d.u)) t.q = mm_quantile(select(d.u, d.u:<.),
                            rows(w)!=1 ? select(w, d.u:<.) : w, t.tau)
    else               t.q = mm_quantile(editmissing(d.u,0), w, t.tau)
    
    // return
    return(t)
}

real scalar mm_aqreg_ymean(struct mm_aqreg_struct scalar t)
{
    if (t.ymean==.z) {
        if (rows(*t.y)==0) t.ymean = .
        else               t.ymean = quadcross(*t.w, *t.y) / mm_aqreg_N(t)
    }
    return(t.ymean)
}

real rowvector mm_aqreg_means(struct mm_aqreg_struct scalar t)
{
    if (t.means==.z) {
        if      (t.k==0)        t.means = J(1, 0, .)
        else if (rows(*t.X)==0) t.means = J(1, t.k, .)
        else                    t.means = quadcross(*t.w, *t.X) / mm_aqreg_N(t)
    }
    return(t.means)
}

real matrix mm_aqreg_b(struct mm_aqreg_struct scalar t)
{
    return((mm_ls_b(t.L) :+ mm_ls_b(t.S) :* J(t.k, 1, t.q)) \ 
           (t.beta0 :+ t.gamma0 :* t.q))
}

real matrix mm_aqreg_xb(struct mm_aqreg_struct scalar t, | real matrix X)
{
    real matrix b
    
    b = mm_aqreg_b(t)
    if (args()==2) {
        if (cols(X)!=t.k) _error(3200)
        if (t.k) return(X * b[|1,1 \ t.k,.|] :+ b[t.k+1,])
        return(J(rows(X), 1, b))
    }
    if (t.k) return(*t.X * b[|1,1 \ t.k,.|] :+ b[t.k+1,])
    return(J(rows(*t.y), 1, b))
}

real matrix mm_aqreg_ue(struct mm_aqreg_struct scalar t)
{
    return(*t.y :- mm_aqreg_xb(t))
}

real matrix mm_aqreg_xbu(struct mm_aqreg_struct scalar t)
{
    return(mm_aqreg_xb(t) + mm_aqreg_u(t))
}

real matrix mm_aqreg_u(struct mm_aqreg_struct scalar t)
{
    return(t.alpha :+ t.delta :* J(rows(t.delta), 1, t.q))
}

real matrix mm_aqreg_e(struct mm_aqreg_struct scalar t)
{
    return(*t.y :- mm_aqreg_xbu(t))
}

real colvector mm_aqreg_levels(struct mm_aqreg_struct scalar t)
{
    return(t.levels)
}

real colvector mm_aqreg_k_levels(struct mm_aqreg_struct scalar t)
{
    return(rows(t.levels))
}

real colvector mm_aqreg_omit(struct mm_aqreg_struct scalar t)
{
    return(mm_ls_omit(t.L) \ 0)
}

real colvector mm_aqreg_k_omit(struct mm_aqreg_struct scalar t)
{
    return(mm_ls_k_omit(t.L))
}

real scalar mm_aqreg_N(struct mm_aqreg_struct scalar t)
{
    if (t.N==.z) {
        t.N = rows(*t.w)==1 ? *t.w * rows(*t.y) : quadsum(*t.w)
    }
    return(t.N)
}

real colvector mm_aqreg_n(struct mm_aqreg_struct scalar t)
{
    return(t.n)
}

real colvector mm_aqreg_beta(struct mm_aqreg_struct scalar t)
{
    return(mm_ls_b(t.L) \ t.beta0)
}

real colvector mm_aqreg_alpha(struct mm_aqreg_struct scalar t)
{
    return(t.alpha)
}

real colvector mm_aqreg_gamma(struct mm_aqreg_struct scalar t)
{
    return(mm_ls_b(t.S) \ t.gamma0)
}

real colvector mm_aqreg_delta(struct mm_aqreg_struct scalar t)
{
    return(t.delta)
}

real rowvector mm_aqreg_tau(struct mm_aqreg_struct scalar t)
{
    return(t.tau)
}

real rowvector mm_aqreg_q(struct mm_aqreg_struct scalar t)
{
    return(t.q)
}

end
