*! version 1.0.2  03mar2021  Ben Jann
version 11
mata:

// M-estimate of location
struct _mm_mls_struct {
    real scalar eff   // not used by mscale
    real scalar bp    // only used by mscale
    real scalar delta // only used by mscale
    real scalar k
    real scalar b
    real scalar b0
    real scalar s     // not used by mscale
    real scalar l     // only used by mscale
    real scalar d
    real scalar tol
    real scalar iter
    real scalar maxiter
    real scalar trace
    real scalar conv
}

struct _mm_mls_struct scalar mm_mloc(real colvector x,
  | real colvector w, real scalar eff, string scalar obj,
    real scalar b, real scalar s, 
    real scalar trace, real scalar tol, real scalar iter)
{
    real scalar b0, i
    real colvector z
    struct _mm_mls_struct scalar S
    pointer(real scalar function) scalar f
    
    if (args()<2) w = 1
    if (hasmissing(w)) _error(3351)
    if (hasmissing(x)) _error(3351)
    S.eff     = (eff<.   ? eff   : 95)
    S.trace   = (trace<. ? trace : 0)
    S.tol     = (tol<.   ? tol   : 1e-10)
    S.maxiter = (iter<.  ? iter  : st_numscalar("c(maxiter)"))
    S.iter    = 0
    S.conv    = 0
    S.b       = .
    if (obj=="" | obj=="huber") {
        f = &mm_huber_w()
        S.k = mm_huber_k(S.eff)
    }
    else if (obj=="biweight") {
        f = &mm_biweight_w()
        S.k = mm_biweight_k(S.eff)
    }
    else {
        printf("{err}'%s' not supported; invalid objective function\n", obj)
        _error(3498)
    }
    if (rows(x)==0) return(S) // nothing to do (b = ., conv = 0)
    S.b = S.b0 = (b<. ? b : mm_median(x, w))
    S.s  = (s<. ? s : mm_median(abs(x :- (b<. ? mm_median(x, w) : S.b)), w)
                      / invnormal(0.75))
    if (S.s<=0) return(S) // scale is 0 (b = starting value, conv = 0)
    // iterate
    z = (x :- S.b) / S.s
    for (i=1; i<=S.maxiter; i++) {
        b0 = S.b
        z = (*f)(z, S.k) :* w
        S.b = mean(x, z)
        S.d = abs(S.b - b0) / S.s
        if (S.trace) printf("{txt}{lalign 16:Iteration %g:}" +
            "absolute difference = %9.0g\n", i, S.d)
        if (S.d <= S.tol) break
        z = (x :- S.b) / S.s
    }
    S.iter = i
    S.conv = (S.iter<=S.maxiter)
    return(S)
}

real scalar    mm_mloc_b(struct _mm_mls_struct scalar S) return(S.b)
real scalar   mm_mloc_b0(struct _mm_mls_struct scalar S) return(S.b0)
real scalar    mm_mloc_s(struct _mm_mls_struct scalar S) return(S.s)
real scalar mm_mloc_conv(struct _mm_mls_struct scalar S) return(S.conv)
real scalar    mm_mloc_d(struct _mm_mls_struct scalar S) return(S.d)
real scalar mm_mloc_iter(struct _mm_mls_struct scalar S) return(S.iter)
real scalar    mm_mloc_k(struct _mm_mls_struct scalar S) return(S.k)
real scalar  mm_mloc_eff(struct _mm_mls_struct scalar S) return(S.eff)

// M-estimate of scale
struct _mm_mls_struct scalar mm_mscale(real colvector x,
  | real colvector w, real scalar bp,
    real scalar b, real scalar l, 
    real scalar trace, real scalar tol, real scalar iter)
{
    real scalar b0, i
    real colvector z
    struct _mm_mls_struct scalar S
    
    if (args()<2) w = 1
    if (hasmissing(w)) _error(3351)
    if (hasmissing(x)) _error(3351)
    S.bp      = (bp<.    ? bp    : 50)
    S.trace   = (trace<. ? trace : 0)
    S.tol     = (tol<.   ? tol   : 1e-10)
    S.maxiter = (iter<.  ? iter  : st_numscalar("c(maxiter)"))
    S.iter    = 0
    S.conv    = 0
    S.b       = .
    S.k       = mm_biweight_k_bp(S.bp)
    if (rows(x)==0) return(S) // nothing to do (b = ., conv = 0)
    S.l = (l<. ? l : mm_median(x, w))
    S.b = S.b0 = (b<. ? b : 
        mm_median(abs(x :- (l<. ? mm_median(x, w) : S.l)), w) / invnormal(0.75))
    if (S.b<=0) return(S) // scale is 0 (b = starting value, conv = 0)
    // iterate
    S.delta = S.bp/100 * S.k^2/6
    z = (x :- S.l)
    for (i=1; i<=S.maxiter; i++) {
        b0 = S.b
        S.b = sqrt(mean(mm_biweight_rho(z/b0, S.k), w) / S.delta) * b0
        S.d = abs(S.b/b0 - 1)
        if (S.trace) printf("{txt}{lalign 16:Iteration %g:}" +
            "absolute relative difference = %9.0g\n", i, S.d)
        if (S.d <= S.tol) break
    }
    S.iter = i
    S.conv = (S.iter<=S.maxiter)
    return(S)
}

real scalar     mm_mscale_b(struct _mm_mls_struct scalar S) return(S.b)
real scalar    mm_mscale_b0(struct _mm_mls_struct scalar S) return(S.b0)
real scalar     mm_mscale_l(struct _mm_mls_struct scalar S) return(S.l)
real scalar  mm_mscale_conv(struct _mm_mls_struct scalar S) return(S.conv)
real scalar     mm_mscale_d(struct _mm_mls_struct scalar S) return(S.d)
real scalar  mm_mscale_iter(struct _mm_mls_struct scalar S) return(S.iter)
real scalar     mm_mscale_k(struct _mm_mls_struct scalar S) return(S.k)
real scalar    mm_mscale_bp(struct _mm_mls_struct scalar S) return(S.bp)
real scalar mm_mscale_delta(struct _mm_mls_struct scalar S) return(S.delta)

// Huber tuning constant for given efficiency
real scalar mm_huber_k(real scalar eff)
{
    if (eff==95) return(1.34499751)
    if (eff==90) return( .98180232)
    if (eff==85) return( .73173882)
    if (eff==80) return( .52942958)
    if (eff<63.7)  _error(3498, "efficiency may not be smaller than 63.7")
    if (eff>99.9) _error(3498, "efficiency may not be larger than 99.9")
    return(round(mm_finvert(eff/100, &mm_huber_eff(), 0.001, 3), 1e-8))
}

// Huber efficiency for given tuning constant
real scalar mm_huber_eff(real scalar k)
{
    return((normal(k)-normal(-k))^2 / 
        (2 * (k^2 * (1 - normal(k)) + normal(k) - 0.5 - k * normalden(k))))
}

// Huber objective functions and weights
real colvector mm_huber_rho(real colvector x, real scalar k)
{
    real colvector y, d
    
    y = abs(x)
    d = y:<=k
    return((0.5 * x:^2):*d :+ (k*y :- 0.5*k^2):*(1 :- d))
}

real colvector mm_huber_psi(real colvector x, real scalar k)
{
    real colvector d
    
    d = abs(x):<=k
    return(x:*d :+ (sign(x)*k):*(1:-d))
}

real colvector mm_huber_phi(real colvector x, real scalar k)
{
    return(abs(x):<= k)
}

real colvector mm_huber_w(real colvector x, real scalar k)
{
    real colvector y
    
    y = abs(x)
    return(editmissing(k :/ y, 0):^(y:>k))
}

// biweight tuning constant for given efficiency
real scalar mm_biweight_k(real scalar eff0)
{
    real scalar k, eff
    
    if (eff0==95) return(4.6850649)
    if (eff0==90) return(3.8826616)
    if (eff0==85) return(3.4436898)
    if (eff0==80) return(3.1369087)
    if (eff0<0.1)  _error(3498, "efficiency may not be smaller than 0.1")
    if (eff0>99.9) _error(3498, "efficiency may not be larger than 99.9")
    eff = eff0/100
    k = 0.8376 + 1.499*eff + 0.7509*sin(1.301*eff^6) + 
        0.04945*eff/sin(3.136*eff) + 0.9212*eff/cos(1.301*eff^6)
    return(round(mm_finvert(eff, &mm_biweight_eff(), k/5, k*1.1), 1e-7))
}

// biweight efficiency for given tuning constant
real scalar mm_biweight_eff(real scalar k)
{   // using Simpson's rule integration
    real scalar    l, u, n, d, phi, psi2
    real colvector x, w
    
    l = 0; u = k; n = 1000; d = (u-l)/n 
    x = rangen(l, u, n+1)
    w = 1 \ colshape((J(n/2,1,4), J(n/2,1,2)), 1)
    w[n+1] = 1
    phi = 2 * (d / 3) * quadcolsum(_mm_biweight_eff_phi(x, k):*w)
    psi2 = 2 * (d / 3) * quadcolsum(_mm_biweight_eff_psi2(x, k):*w)
    return(phi^2 / psi2)
}

real matrix _mm_biweight_eff_phi(real matrix x, real scalar k)
{
    real matrix x2

    x2 = (x / k):^2
    return(normalden(x) :* ((1 :- x2) :* (1 :- 5*x2)) :* (x2:<=1))
}

real matrix _mm_biweight_eff_psi2(real matrix x, real scalar k)
{
    real matrix x2

    x2 = (x / k):^2
    return(normalden(x) :* ((x :* (1 :- x2):^2) :* (x2:<=1)):^2)
}

// biweight tuning constant for given breakdown point
real scalar mm_biweight_k_bp(real scalar bp)
{
    if (bp==50) return(1.547645)
    if (bp<1)  _error(3498, "bp may not be smaller than 1")
    if (bp>50) _error(3498, "bp may not be larger than 50")
    return(round(mm_finvert(bp/100, &mm_biweight_bp(), 1.5, 18), 1e-7))
}

// biweight breakdown point for given tuning constant
real scalar mm_biweight_bp(real scalar k)
{   // using Simpson's rule integration
    real scalar    l, u, n, d
    real colvector x, w

    l = 0; u = k; n = 1000; d = (u-l)/n
    x = rangen(l, u, n+1)
    w = 1 \ colshape((J(n/2,1,4), J(n/2,1,2)), 1)
    w[n+1] = 1
    return(2 * (normal(-k) +  d/3 *
        quadcolsum(normalden(x):*mm_biweight_rho(x, k):*w) / (k^2/6)))
}

// biweight objective functions and weights
real colvector mm_biweight_rho(real colvector x, real scalar k)
{
    real colvector x2

    x2 = (x / k):^2
    return(k^2/6 * (1 :- (1 :- x2):^3):^(x2:<=1))
}

real colvector mm_biweight_psi(real colvector x, real scalar k)
{
    real colvector x2

    x2 = (x / k):^2
    return((x :* (1 :- x2):^2) :* (x2:<=1))
}

real colvector mm_biweight_phi(real colvector x, real scalar k)
{
    real colvector x2

    x2 = (x / k):^2
    return(((1 :- x2) :* (1 :- 5*x2)) :* (x2:<=1))
}
real colvector mm_biweight_w(real colvector x, real scalar k)
{
    real colvector x2

    x2 = (x / k):^2
    return(((1 :- x2):^2) :* (x2:<=1))
}

end
