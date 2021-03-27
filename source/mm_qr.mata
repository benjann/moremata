*! version 1.0.1  27mar2021  Ben Jann
version 11
mata:

class mm_qr
{
    // constructors
    private:
        void           new()
        void           init() // default settings
        void           clear(), clear1(), clear2(), clear2b(), clear3()

    // setup
    public:
        void           data()
        transmorphic   qd()
        transmorphic   demean()
        transmorphic   collin()
        transmorphic   p()
        transmorphic   b_init()
        transmorphic   tol()
        transmorphic   maxiter()
        transmorphic   beta()
        transmorphic   log()
    
    private:
        pointer scalar y
        pointer scalar X
        pointer scalar w
        real scalar    n, N
        real scalar    cons
        real scalar    k, kadj
        real scalar    K, Kadj
        real scalar    p
        real colvector b_init, b_ls
        real scalar    b_init_user
        real scalar    maxiter
        real scalar    tol
        real scalar    beta
        real scalar    qd
        real scalar    demean
        real scalar    collin
        real colvector omit
        real rowvector xindx
        real scalar    log
    
    // results
    public:
        real colvector b()
        real colvector xb()
        real scalar    gap()
        real scalar    sdev()
        real scalar    iter()
        real scalar    converged()
        real scalar    n() // number of observations
        real scalar    N() // sum of weights
        real scalar    cons()
        real scalar    k() // number of predictors
        real scalar    K() // number of coefficients
        real colvector omit()
        real scalar    k_omit()
        real scalar    ymean()
        real rowvector means()
    
    private:
        real colvector b
        real scalar    conv
        real scalar    iter
        real scalar    gap
        real scalar    sdev
        real scalar    ymean
        real rowvector means
        
    // functions
    private:
        void           err_nodata()
        void           printflush()
        void           printlog()
        real colvector _xb()
        real scalar    _sdev()
        void           set_collin()
        transmorphic   lsfit()
        void           _b_init()
        real colvector rmomit(), addomit()
        real colvector meanadj()
        void           fit(), fnb()
        pointer scalar get_X()
        real colvector unzero()
        real matrix    cross()
        real scalar    minselect()
}

void mm_qr::new()
{
    init()
    clear()
}

void mm_qr::init()
{
    p         = .5
    qd        = 1
    demean    = 1
    collin    = 1
    tol       = 1e-10
    maxiter   = st_numscalar("c(maxiter)")
    beta      = 0.99995
    log       = 0
}

void mm_qr::clear() // if new() or data()
{
    y = X = w = NULL
    n = N = cons = k = K = ymean = means = .z
    b_init_user = 0 // b_init can only be set after data has been set
    clear1()
}

void mm_qr::clear1() // if qd() or demean()
{
    b_ls = .z
    if (!b_init_user) b_init = .z
    clear2()
}

void mm_qr::clear2() // if collin()
{
    kadj = Kadj = omit = xindx = .z
    clear3()
}

void mm_qr::clear2b() // if p()
{
    if (!b_init_user & cons) b_init = .z
    clear3()
}

void mm_qr::clear3() // if b_init(), tol(), maxiter(), beta()
{
    b = conv = iter = gap = sdev = .z
}

void mm_qr::err_nodata()
{
    if (y==NULL) {
        display("{err}data not set")
        _error(3498)
    }
}

void mm_qr::printflush(string scalar s)
{
    printf(s)
    displayflush()
}

void mm_qr::printlog(real scalar d, real scalar iter)
{
    printf("{txt}")
    printf("{txt}Iteration %g:", iter)
    printf("{col 16}mreldif() in b = {res}%11.0g;", d)
    printf("{txt}  duality gap = {res}%11.0g\n", gap)
    displayflush()
}

void mm_qr::data(real colvector y0, | real matrix X0, real colvector w0,
    real colvector cons0)
{
    clear()
    // do error checks before storing anything
    if (!(X0==. | X0==J(0,0,.))) {
        if (rows(X0)!=rows(y0)) _error(3200)
    }
    if (rows(w0)) {
        if (rows(w0)!=rows(y0) & rows(w0)!=1) _error(3200)
    }
    // now start storing
    y = &y0
    n = rows(y0)
    if (X0==. | X0==J(0,0,.)) X = &J(n, 0, 1)
    else                      X = &X0
    if (rows(w0)) {
        if (rows(w0)==1) N = n * w0
        else             N = quadsum(w0)
        w = &w0
    }
    else {
        N = n
        w = &1
    }
    cons = (cons0!=0)
    k = cols(*X)
    K = k + cons
}

transmorphic mm_qr::qd(| real scalar qd0)
{
    if (args()==0) return(qd)
    if (qd==(qd0!=0)) return // no change
    qd = (qd0!=0)
    clear1()
}

transmorphic mm_qr::demean(| real scalar demean0)
{
    if (args()==0) return(demean)
    if (demean==(demean0!=0)) return // no change
    demean = (demean0!=0)
    clear1()
}

transmorphic mm_qr::collin(| real scalar collin0)
{
    if (args()==0) return(collin)
    if (collin==(collin0!=0)) return // no change
    collin = (collin0!=0)
    clear2()
}

transmorphic mm_qr::p(| real scalar p0)
{
    if (args()==0) return(p)
    if (p0<=0 | p0>=1) _error(3300)
    if (p==p0) return // no change
    p = p0
    clear2b()
}

transmorphic mm_qr::b_init(| real colvector b_init0)
{
    if (args()==0) {
        if (b_init==.z) {
            err_nodata()
            _b_init()
        }
        return(b_init)
    }
    if (b_init0==.z) {  // use b_init(.z) to clear starting values
        b_init_user = 0
    }
    else {
        err_nodata()
        if (rows(b_init0)!=K) _error(3200) // wrong number of parameters
        b_init = b_init0
        b_init_user = 1
    }
    clear3()
}

transmorphic mm_qr::tol(| real scalar tol0)
{
    if (args()==0) return(tol)
    if (tol0<=0) _error(3300)
    if (tol==tol0) return // no change
    tol = tol0
    clear3()
}

transmorphic mm_qr::maxiter(| real scalar maxiter0)
{
    if (args()==0) return(maxiter)
    if (maxiter0<0) _error(3300)
    if (maxiter==maxiter0) return // no change
    maxiter = maxiter0
    clear3()
}

transmorphic mm_qr::beta(| real scalar beta0)
{
    if (args()==0) return(beta)
    if (beta0<0 | beta0>1) _error(3300)
    if (beta==beta0) return // no change
    beta = beta0
    clear3()
}

transmorphic mm_qr::log(| real scalar log0)
{
    if (args()==0) return(log)
    log = log0
}

real scalar mm_qr::n() return(n)

real scalar mm_qr::N() return(N)

real scalar mm_qr::cons() return(cons)

real scalar mm_qr::k() return(k)

real scalar mm_qr::K() return(K)

real colvector mm_qr::omit()
{
    if (omit==.z) {
        if (y!=NULL) set_collin()
    }
    return(omit)
}

real scalar mm_qr::k_omit()
{
    if (kadj==.z) {
        if (y!=NULL) set_collin()
    }
    return(K-Kadj)
}

real scalar mm_qr::ymean()
{
    if (ymean==.z) {
        if (y!=NULL) ymean = quadcross(*w, *y) / N
    }
    return(ymean)
}

real rowvector mm_qr::means()
{
    if (means==.z) {
        if (y!=NULL) means = quadcross(*w, *X) / N
    }
    return(means)
}

real colvector mm_qr::b()
{
    if (b==.z) {
        err_nodata()
        fit()
    }
    return(b)
}

real colvector mm_qr::xb(| real matrix X0)
{
    if (b==.z) (void) b()
    if (args()==1) {
        if (cols(X0)!=k) _error(3200)
        return(_xb(X0, b))
    }
    return(_xb(*X, b))
}

real colvector mm_qr::_xb(real matrix X, real colvector b)
{
    real scalar k
    
    k = rows(b) - cons
    if (k==0) {
        if (cons) return(J(rows(*y), 1, b))
        return(J(rows(*y), 1, .)) // model has no parameters
    }
    if (cons) return(X * b[|1\k|] :+ b[k+1])
    return(X * b)
}

real scalar mm_qr::gap()
{
    if (b==.z) (void) b()
    return(gap)
}

real scalar mm_qr::sdev()
{
    if (b==.z) (void) b()
    if (sdev==.z) sdev = _sdev(b)
    return(sdev)
}

real scalar mm_qr::_sdev(real colvector b)
{
    real colvector e
    
    e = *y - _xb(*X, b)
    return(cross(*w, (p :- (e:<0)) :* e))
}

real scalar mm_qr::iter()
{
    if (b==.z) (void) b()
    return(iter)
}

real scalar mm_qr::converged()
{
    if (b==.z) (void) b()
    return(conv)
}

void mm_qr::fit()
{
    real colvector b
    
    // identify collinear variables and exit if nothing to do
    set_collin()
    iter = conv = 0
    if (Kadj==0) { // model has no (non-omitted) parameters
        this.b = J(K,1,.)
        return
    }
    if (n==0) { // no observations
        this.b = J(K,1,.)
        return
    }
    // starting values
    if (b_init==.z) _b_init()
    b = -rmomit(meanadj(b_init, -1))
    
    // interior point algorithm (b will be replaced by solution)
    fnb(b)
    
    // rescale coefficients
    this.b = meanadj(addomit(-b), 1)
}

void mm_qr::set_collin()
{
    real scalar  nomit
    transmorphic t
    
    if (kadj!=.z) return    // already applied
    if (collin==0 | k==0) {
        kadj  = k
        Kadj  = K
        omit  = J(K, 1, 0)
        xindx = J(1, 0, .)
        return
    }
    t = lsfit()
    omit = mm_ls_omit(t)
    if (cons) omit[K] = 0 // not necessary, I believe
    nomit = sum(omit)
    kadj = k - nomit
    Kadj = K - nomit
    if (!nomit) xindx = J(1, 0, .)
    else        xindx = select(1..k, !omit[|1\k|]')
}

transmorphic mm_qr::lsfit()
{
    transmorphic t
    
    t = mm_ls(*y, *X, *w, cons, qd, demean)
    // hold on to coefficients for b_init
    b_ls = mm_ls_b(t)
    // hold on to means for demeaning
    if (cons & demean) { 
        ymean = mm_ls_ymean(t)
        means = mm_ls_means(t)
    }
    return(t)
}

void mm_qr::_b_init()
{
    if (K==0) { // model has no parameters
        b_init = J(K,1,.)
        return
    }
    if (n==0) { // no observations
        b_init = J(K,1,.)
        return
    }
    if (b_ls==.z) (void) lsfit()
    b_init = b_ls
    if (cons) b_init[K] = b_init[K] + mm_quantile(*y - _xb(*X, b_ls), *w, 1-p)
}

real colvector mm_qr::rmomit(real colvector b)
{
    if (length(xindx)==0) return(b)
    return(b[xindx] \ (cons ? b[K] : J(0,1,.)))
}

real colvector mm_qr::addomit(real colvector b0)
{
    real colvector b
    
    if (length(xindx)==0) return(b0)
    b = J(K,1,0)
    b[xindx] = b0[|1\kadj|]
    if (cons) b[K] = b0[Kadj]
    return(b)
}

real colvector mm_qr::meanadj(real colvector b0, real scalar sign)
{
    real colvector b
    
    if (!demean) return(b0)
    if (!cons)   return(b0)
    b = b0
    b[K] = b[K] + ymean() * sign
    if (k) b[K] = b[K] - (means()*b[|1\k|]) * sign
    return(b)
}

// loose translation of rqfnb.f from "quantreg" package version 5.85 for R
// - invsym() is used inside the algorithm instead of a solver such as
//   lusolve() because invsym() appears more robust against failure due to
//   numerical problems and because accuracy does not seem to be a major
//   concern within the iterations
void mm_qr::fnb(real colvector y)
{
    real scalar    fp, fd, mu, g
    real colvector dy, y0
    real colvector W, x, s, z, w, dx, ds, dz, dw, dxdz, dsdw, dr, d, dW
    real matrix    iada
    pointer scalar a
    
    // obtain data; a copy of X will only be made if necessary; for simplicity,
    // weights are normalized to the number of observations
    s = -(*this.y)
    if (cons & demean) s = s :+ ymean()
    a = get_X()
    if (rows(*this.w)!=1) W = (*this.w) * (n/N)
    else                  W = 1
    
    // compute residuals to be used as input for the algorithm; if necessary,
    // the vector of residuals will be shifted up or down by a small amount so
    // that values equal to or very close to 0 (which would make the algorithm 
    // fail) are avoided; in this case, the initial coefficients (y) will be
    // adjusted so that they are consistent with the shifted residuals (to
    // avoid bias due to an implicit change in the dependent variable); the
    // procedure may fail in which case residuals close to zero will moved by a
    // fixed amount without updating the initial coefficients (this may lead to
    // a small bias in the estimates)
    s = unzero(s, *a, y)

    // algorithm
    z = s :* (s :> 0)
    w = z - s
    x = J(n, 1, 1-p)
    s = J(n, 1, p)
    if (maxiter<=0) {
        // make sure gap is not empty even if maxiter=0
        gap = cross(z, W, x) + cross(w, W, s)
    }
    while (iter < maxiter) {
        iter = iter + 1
        d    = 1 :/ (z :/ x + w :/ s)
        dW   = W :* (d / max(d)) // set max(d)=1 for numerical stability
        dr   = z - w
        iada = invsym(cross(*a,cons, dW, *a,cons))
        dy   = iada * cross(*a,cons, dW, dr,0)
        dx   = d :* (_xb(*a, dy) - dr)
        ds   = -dx
        dz   = -z :* (1 :+ dx :/ x)
        dw   = -w :* (1 :+ ds :/ s)
        fp   = minselect(x, dx, s, ds)
        fd   = minselect(w, dw, z, dz)
        if (min((fp, fd)) < 1) {
            mu   = cross(z, W, x) + cross(w, W, s)
            g    = cross(z + fd * dz, W, x + fp * dx) + 
                   cross(w + fd * dw, W, s + fp * ds)
            mu   = mu * (g / mu)^3 / (2 * n)
            dxdz = dx :* dz
            dsdw = ds :* dw
            dr   = dr + mu:/s - mu:/x + dxdz:/x - dsdw:/s
            dy   = iada * cross(*a,cons, dW, dr,0)
            dx   = d :* (_xb(*a, dy) - dr)
            ds   = -dx
            dz   = (mu :- z :* dx :- dxdz) :/ x - z
            dw   = (mu :- w :* ds :- dsdw) :/ s - w
            fp   = minselect(x, dx, s, ds)
            fd   = minselect(w, dw, z, dz)
        }
        x    = x + fp * dx
        s    = s + fp * ds
        w    = w + fd * dw
        z    = z + fd * dz
        y0   = y
        y    = y + fd * dy
        gap = cross(z, W, x) + cross(w, W, s)
        if (log) {
            if (log==2) printflush(".")
            else        printlog(mreldif(y, y0), iter)
        }
        if (gap<tol) {
            conv = 1
            break
        }
    }
    if (log==2 & iter) printflush("\n")
}

pointer scalar mm_qr::get_X()
{
    real matrix x
    
    if (!(cons & demean) & !length(xindx)) return(X)
    x = *X
    if (cons & demean) x = x :- means()
    if (length(xindx)) x = x[,xindx]
    return(&x)
}

real colvector mm_qr::unzero(real colvector y, real matrix X, real colvector b)
{
    real scalar    eps, omax
    real scalar    lo, up, offset
    real colvector r, rlo, rup, b1
    
    // settings
    eps  = 1e-06
    omax = 0.001
    // obtain residuals from initial coefficients and check whether any of them
    // are equal to (or very close to) zero
    r   = y - _xb(X, b)
    rlo = (r :< -eps)
    rup = (r :>  eps)
    if (all(rlo :+ rup)) {
        // no residuals within [-eps,eps]; ok to proceed
        return(r)
    }
    // compute reasonable offset; we use half the distance to closest residual
    // outside [-eps,eps]
    lo = max(select(r, rlo)) + eps
    up = min(select(r, rup)) - eps
    if (up<-lo) offset = up/2
    else        offset = lo/2
    if (offset>=.) offset = omax // can happen in case of perfect fit
    // modify residuals and initial coefficients (both have to be modified so
    // that the implied values of the dependent variable do not change)
    if (cons) {
        // if model has a constant: simply shift the residuals and then adjust
        // the constant
        r = r :+ offset
        if (abs(offset)>2*eps) {
            // shift cannot fail if abs(offset)>2*eps
            b[cols(X)+1] = b[cols(X)+1] - offset
            return(r)
        }
        if (all((r :< -eps) :| (r :>  eps))) {
            // shift was successful despite abs(offset)<=2*eps
            b[cols(X)+1] = b[cols(X)+1] - offset
            return(r)
        }
    }
    else {
        // if model has no constant:
        // (1) regress X on offset
        // (2) update initial coefficients using coefficients from (1)
        // (3) recompute residuals based on updated coefficients
        b1 = mm_lsfit(J(n,1,offset), X, *w, 0) // (weights not strictly needed)
        b1 = b - b1
        r = y - _xb(X, b1)
        if (all((r :< -eps) :| (r :>  eps))) {
            // shift was successful
            b = b1
            return(r)
        }
    }
    // fallback rule: only move critical residuals away from zero; do not update
    // coefficients; this may lead to somewhat biased estimates; a better
    // approach would be to return to the above approach and keep on searching
    // for suitable starting values
    r = y - _xb(X, b)
    offset = min((omax, max((eps,abs(offset))))) // offset in [eps, omax]
    return(r - offset * (r:>-eps :& r:<0) + offset * (r:>=0 :& r:<eps))
}

real matrix mm_qr::cross(real matrix a, real matrix b, | real matrix c,
    real matrix d, real matrix e)
{
    if (args()==2) {
        if (qd) return(quadcross(a, b))
        return(::cross(a, b))
    }
    if (args()==3) {
        if (qd) return(quadcross(a, b, c))
        return(::cross(a, b, c))
    }
    if (args()==4) {
        if (qd) return(quadcross(a, b, c, d))
        return(::cross(a, b, c, d))
    }
    if (qd) return(quadcross(a, b, c, d, e))
    return(::cross(a, b, c, d, e))
}

real scalar mm_qr::minselect(real colvector x, real colvector dx, 
    real colvector s, real colvector ds)
{
    real colvector min, p
    
    min = J(3,1,.)
    p = select(1::n, dx:<0)
    if (length(p)) min[1] = min(-x[p] :/ dx[p]) * beta
    p = select(1::n, ds:<0)
    if (length(p)) min[2] = min(-s[p] :/ ds[p]) * beta
    min[3] = 1
    return(min(min))
}

// - wrapper for quick QR fit
real colvector mm_qrfit(
    real colvector      y,
    | real matrix       X, 
      real colvector    w, 
      real scalar       p,
      real scalar       cons,
      real colvector    b_init,
      real scalar       relax)
{
    class mm_qr scalar S
    
    S.data(y, X, w, cons)
    if (p!=.) S.p(p)
    if (args()>=6 & b_init!=.) S.b_init(b_init)
    if (args()<7) relax = 0
    if (!S.converged() & !relax) _error(3360)
    return(S.b())
}

end
