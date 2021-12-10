*! version 1.0.1  18sep2021  Ben Jann
version 9.2
mata:

// Main functions

real matrix mm_aqregfit(real colvector y, real colvector id, | real matrix X,
    real colvector w, real vector tau, real scalar sort, real scalar qd)
{
    if (args()<3) X = .
    if (args()<4) w = 1
    if (args()<5) tau = 0.5
    return(mm_aqreg_b(mm_aqreg(y, id, X, w, tau, sort, qd)))
}

struct mm_aqreg_struct {
    pointer(real colvector) scalar y, w
    pointer(real matrix) scalar    X
    real scalar                    k, N, ymean
    real rowvector                 means
    real rowvector                 tau, q
    real scalar                    beta0, gamma0
    real colvector                 alpha, delta
    struct mm_ls_struct scalar     L, S
    pointer(struct mm_areg_struct_g scalar) scalar g
}

struct mm_aqreg_struct scalar mm_aqreg(
    real colvector   y,
    real colvector   id,
    | real matrix    X,
      real colvector w,
      real vector    tau,
      real scalar    sort,
      real scalar    qd)
{
    if (args()<3) X = .
    if (args()<4) w = 1
    if (args()<5) tau = 0.5
    if (rows(id)!=rows(y)) _error(3200)
    return(_mm_aqreg(_mm_areg_g(id, sort), y, X, w, tau, qd))
}

struct mm_aqreg_struct scalar _mm_aqreg(
    struct mm_areg_struct_g scalar g,
    real colvector                 y,
    | real matrix                  X,
      real colvector               w,
      real vector                  tau,
      real scalar                  qd)
{
    real colvector ym, yd, r, rm, rd, e, u
    real matrix    Xm, Xd
    struct mm_aqreg_struct scalar t
    
    // setup
    if (args()<3) X = .
    if (args()<4) w = 1
    if (args()<5) tau = 0.5
    t.N = t.ymean = t.means = .z
    t.g = &g
    t.y = &y
    if (X==.) t.X = &(J(rows(y), 0, .))
    else      t.X = &X
    t.w = &w
    t.k = cols(*t.X)
    
    // quantiles
    if (rows(tau)>1) t.tau = tau'
    else             t.tau = tau
    if (length(tau)==0) _error(3200)
    if (any(tau:<=0 :| tau:>=1)) _error(3300, "tau must be between 0 and 1")
    
    // transform data
    ym = y; _mm_areg_gmean(g, ym, w, qd); yd = y - ym
    if (t.k) {; Xm = *t.X; _mm_areg_gmean(g, Xm, w, qd); Xd = *t.X - Xm; }
    else        Xd = Xm = *t.X
    
    // location fit (beta0: global constant; alpha: fixed effects)
    t.L = mm_ls(yd, Xd, w, 0, qd)
    t.L.y = t.L.X = t.L.w = NULL // release pointers to data (save memory)
    if (t.k) t.beta0 = mm_aqreg_ymean(t) - mm_ls_xb(t.L, mm_aqreg_means(t))
    else     t.beta0 = mm_aqreg_ymean(t)
    if (t.k) t.alpha = ym - mm_ls_xb(t.L, Xm)
    else     t.alpha = ym
    t.alpha = t.alpha :- t.beta0
    if (t.k) e = yd - mm_ls_xb(t.L, Xd)
    else     e = yd
    
    // scale fit (gamma0: global constant; delta: fixed effects)
    r = 2 * e :* ((e:>=0) :- mean(e:>=0, w))
    rm = r; _mm_areg_gmean(g, rm, w, qd); rd = r - rm
    t.S = mm_ls(rd, Xd, w, 0, qd)
    t.S.y = t.S.X = t.S.w = NULL // release pointers to data (save memory)
    if (t.k) t.gamma0 = mean(r, w) - mm_ls_xb(t.S, mm_aqreg_means(t))
    else     t.gamma0 = mean(r, w)
    if (t.k) t.delta = rm - mm_ls_xb(t.S, Xm)
    else     t.delta = rm
    if (t.k) u = e :/ (mm_ls_xb(t.S, *t.X) + t.delta)
    else     u = e :/ t.delta
    t.delta = t.delta :- t.gamma0
    
    // quantiles
    if (missing(u)) t.q = mm_quantile(select(u, u:<.),
                            rows(w)!=1 ? select(w, u:<.) : w, t.tau)
    else            t.q = mm_quantile(editmissing(u,0), w, t.tau)
    
    // return
    return(t)
}

// Results

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

real colvector mm_aqreg_levels(struct mm_aqreg_struct scalar t)
{
    return(t.g->levels)
}

real colvector mm_aqreg_k_levels(struct mm_aqreg_struct scalar t)
{
    return(rows(t.g->idx))
}

real colvector mm_aqreg_n(struct mm_aqreg_struct scalar t)
{
    return(t.g->n)
}

real rowvector mm_aqreg_tau(struct mm_aqreg_struct scalar t)
{
    return(t.tau)
}

real rowvector mm_aqreg_q(struct mm_aqreg_struct scalar t)
{
    return(t.q)
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

end
