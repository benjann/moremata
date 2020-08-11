*! version 1.0.0  10aug2020  Ben Jann
version 9.2
mata:

real matrix mm_ddens(real colvector X, | real colvector w, real vector minmax0,
    real scalar n0, real scalar h0, real scalar qui) // h0 will be replaced
{
    real scalar    n, N, h, lo, up
    real rowvector minmax
    real colvector AT, W, a
    
    // defaults
    if (args()<2) w  = 1
    if (n0>=.) n0 = 2^14
    if (h0<=0) _error(3300)
    if (args()<6) qui = 0
    
    // check input
    if (missing(X) | missing(w)) _error(3351)
    if (any(w:<0)) {
        display("{err}{it:w} must not be negative")
        _error(3300)
    }
    if (n0<2) _error(3300)
    
    // evaluation grid
    n = 2^ceil(ln(n0)/ln(2)) // round up to next power of 2
    if (length(minmax0)>=1) lo = minmax0[1]
    if (length(minmax0)>=2) up = minmax0[2]
    minmax = minmax(X)
    if (lo>=.) lo = minmax[1] - (minmax[2]-minmax[1])/10 // min - 10% of range
    if (up>=.) up = minmax[2] + (minmax[2]-minmax[1])/10 // min + 10% of range
    if (lo>minmax[1]) {
        display("{err}data smaller than lower bound not allowed")
        _error(3300)
    }
    if (up<minmax[2]) {
        display("{err}data larger than upper bound not allowed")
        _error(3300)
    }
    // bin data
    N  = mm_nobs(X, w)
    AT = rangen(lo, up, n)
    W  = mm_exactbin(X, w, rangen(lo, up, n+1)) / N
        // need to use exact binning because linear binning would introduce
        // some (non-vanishing) bias at the boundaries (doubling the first and
        // last grid count does not seem to help); a consequence of exact 
        // binning is that the density estimate will be slightly shifted/stretched
        // to the left; this error can be substantial if the grid size is small,
        // but it vanishes with increasing grid size
    // obtain discrete cosine transform of binned data
    a = Re( (1 \ 2 * exp(1i * (1::n-1) * pi() / (2*n)))
         :* fft(W[mm_seq(1,n-1,2)] \ W[mm_seq(n,2,2)]) )
    // compute bandwidth (unless provided by user)
    if (h0<.) h = (h0 / (up-lo))^2
    else {
        h = _mm_ddens_h(AT, W, N, (1::n-1):^2, (a[2::n]/2):^2, qui)
        h0 = sqrt(h) * (up-lo) // return optimal bandwidth
    }
    // smooth discrete cosine transform using bandwidth and back-transform
    a = a :* exp(-(0::n-1):^2 * pi()^2 * h/2)
    a = Re(invfft((n * exp(-1i * (0::n-1) * pi() / (2*n))) :* a)) / (up-lo)
    a[mm_seq(1, n, 2) \ mm_seq(2, n, 2)] = a[1::n/2 \ n::n/2+1] // reorder
    // return density estimate and grid
    return((a, AT))
}

real scalar _mm_ddens_h(real colvector AT, real colvector W, real scalar N, 
    real colvector I, real colvector a2, real scalar qui)
{
    real scalar s, hmin, h_os, ax, bx, rc, h
    
    // compute oversmoothed bandwidth
    s = _mm_iqrange(AT, W) / 1.349
    if (s<=0) s = sqrt(variance(AT, W*N))
    else      s = min((s, sqrt(variance(AT, W*N))))
    h_os = (s * (243/(35*N))^.2 * mm_kdel0_gaussian() / (AT[rows(AT)]-AT[1]))^2
    // minimum bandwidth given grid size
    hmin = (.5/(rows(AT)-1) * mm_kdel0_gaussian()/mm_kdel0_rectangle())^2
    // find optimal h
    bx = h_os
    ax = max((hmin, bx/2))
    rc = mm_root(h=., &_mm_ddens_fixpnt(), ax, bx, 0, 100, N, I, a2)
    if (rc==2) { // move down
        bx = ax
        while (1) {
            if (bx<=hmin) break // cannot go below hmin
            ax = max((hmin, bx/2))
            rc = mm_root(h, &_mm_ddens_fixpnt(), ax, bx, 0, 100, N, I, a2)
            if (rc==2) bx = ax // continue moving down
            else break // this also stops if there is a change in 
                      //  direction (rc==3) (in this case: h = current bx)
        }
    }
    else if (rc==3) { // move up
        ax = bx; bx = ax*1.5
        while (1) {
            rc = mm_root(h, &_mm_ddens_fixpnt(), ax, bx, 0, 100, N, I, a2)
            if (rc==3) {; ax = bx; bx = ax*1.5; } // continue moving up
            else break // this also stops if there is a change in 
                       //  direction (rc==2) (in this case: h = current ax)
        }
    }
    if (h<=hmin) { // solution is smaller than hmin
        ax = h_os; bx = ax*1.5
        rc = mm_root(h, &_mm_ddens_fixpnt(), ax, bx, 0, 100, N, I, a2)
        if (rc==2) h = . // moving up does not help
        else if (rc==3) { // continue moving up
            ax = bx; bx = ax*1.5
            while (1) {
                rc = mm_root(h, &_mm_ddens_fixpnt(), ax, bx, 0, 100, N, I, a2)
                if (rc==3) {; ax = bx; bx = ax*1.5; } // continue moving up
                else break // this also stops if there is a change in 
                           //  direction (rc==2) (in this case: h = current ax)
            }
        }
    }
    // return result
    if (h>=.) {
        if (qui==0) display("{txt}(bandwidth estimation failed; using oversmoothed rule)")
        return(h_os)
    }
    return(h)
}

real scalar _mm_ddens_fixpnt(real scalar h, real scalar N, real colvector I, 
    real colvector a2)
{
    real scalar    l, s, K0, c
    real colvector f, t
    
    l = 7
    f = 2 * pi()^(2*l) * sum(I:^l :* a2 :* exp(-I * pi()^2 * h))
    for (s=l-1; s>=2; s--) {
        K0 = mm_prod(mm_seq(1, 2*s-1, 2)) / sqrt(2*pi())
        c  = (1 + (1/2)^(s + 1/2)) / 3
        t  = (2 * c * K0/N/f):^(2/(3 + 2*s))
        f  = 2 * pi()^(2*s) * sum(I:^s :* a2 :* exp(-I * pi()^2 * t))
    }
    return((2 * N * sqrt(pi()) * f)^(-2/5) - h)
}

end

