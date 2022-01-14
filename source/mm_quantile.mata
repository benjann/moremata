*! version 2.0.6  14jan2022  Ben Jann
version 9.2
mata:

real matrix mm_quantile(real matrix X, | real colvector w,
    real matrix P, real scalar d, real scalar fw, real scalar wd)
{
    real colvector o
    pointer (real matrix) scalar XX, ww
    
    if (args()<2) w = 1
    if (args()<3) P = (0, .25, .50, .75, 1)'
    if (args()<4) d = 2
    if (args()<5) fw = 0
    if (!anyof(0::11, d)) {
        display("{err}{it:def} must be an integer in [0,11]")
        _error(3300)
    }
    if (missing(X) | missing(w)) _error(3351)
    if (any(w:<0)) {
        display("{err}{it:w} must not be negative")
        _error(3300)
    }
    // handle weights
    XX = &X; ww = &w
    if (rows(w)!=1) {
        if (rows(w)!=rows(X)) _error(3200)
        // drop zero frequency observations
        if (anyof(w,0)) {
            o = select(1::rows(w), w)
            if (cols(o)==0) o = J(0,1,.) // select() may return 0x0
            XX = &(X[o,])
            ww = &(w[o])
        }
        // weights can be ignored if constant and d < 3
        if (d < 3 & rows(*ww)) { // (*ww may be void)
            if (allof(*ww, (*ww)[1])) ww = &1
        }
        // weights can be ignored if w = 1 for all observations
        else if (allof(*ww,1)) ww = &1
    }
    else if (w==0) {
        XX = &(J(0,cols(X),.))
        ww = &(J(0,1,.))
    }
    // compute quantiles
    if (cols(X)==1 & cols(P)!=1 & rows(P)==1)
        return(_mm_quantile_sort(*XX, *ww, P', d, fw, wd)')
    if (cols(X)!=1 & cols(P)!=1 & cols(X)!=cols(P)) _error(3200)
    return(_mm_quantile_sort(*XX, *ww, P, d, fw, wd))
}

real matrix _mm_quantile_sort(real matrix X, real colvector w,
    real matrix P, real scalar d, real scalar fw, real scalar wd)
{
    real scalar    i, c, c1, c2
    real colvector p, sX, sw, pP, sP
    real matrix    Q

    c1 = cols(X); c2 = cols(P)
    c = max((c1,c2))
    Q = J(rows(P), c, .)
    if (c1==c2) {
        if (w==1) {
            for (i=c; i; i--) {
                Q[,i] = _mm_quantile_d(sort(X[,i],1), editmissing(P[,i],1), d, wd)
            }
            return(Q)
        }
        if (rows(w)==1) {
            for (i=c; i; i--) {
                pP = order(P[,i],1); sP = P[pP,i]
                Q[pP,i] = _mm_quantile_w(sort(X[,i],1), w, sP, d, fw, wd)
            }
            return(Q)
        }
        for (i=c; i; i--) {
            p = order((X[,i],w),(1,2)); sX = X[p,i]; sw = w[p]
            pP = order(P[,i],1); sP = P[pP,i]
            Q[pP,i] = _mm_quantile_w(sX, sw, sP, d, fw, wd)
        }
        return(Q)
    }
    if (c1==1) {
        if (w==1) {
            sX = sort(X,1)
            for (i=c; i; i--) {
                Q[,i] = _mm_quantile_d(sX, editmissing(P[,i],1), d, wd)
            }
            return(Q)
        }
        if (rows(w)==1) {
            sX = sort(X,1)
            for (i=c; i; i--) {
                pP = order(P[,i],1); sP = P[pP,i]
                Q[pP,i] = _mm_quantile_w(sX, w, sP, d, fw, wd)
            }
            return(Q)
        }
        p = order((X,w),(1,2)); sX = X[p]; sw = w[p]
        for (i=c; i; i--) {
            pP = order(P[,i],1); sP = P[pP,i]
            Q[pP,i] = _mm_quantile_w(sX, sw, sP, d, fw, wd)
        }
        return(Q)
    }
    if (c2==1) {
        if (w==1) {
            sP = editmissing(P,1)
            for (i=c; i; i--) {
                Q[,i] = _mm_quantile_d(sort(X[,i],1), sP, d, wd)
            }
            return(Q)
        }
        if (rows(w)==1) {
            pP = order(P,1); sP = P[pP]
            for (i=c; i; i--) {
                Q[pP,i] = _mm_quantile_w(sort(X[,i],1), w, sP, d, fw, wd)
            }
            return(Q)
        }
        pP = order(P,1); sP = P[pP]
        for (i=c; i; i--) {
            p = order((X[,i],w),(1,2)); sX = X[p,i]; sw = w[p]
            Q[pP,i] = _mm_quantile_w(sX, sw, sP, d, fw, wd)
        }
        return(Q)
    }
    _error(3200)
}

real matrix _mm_quantile(real colvector X, | real colvector w, 
    real matrix P, real scalar d, real scalar fw, real scalar wd)
{   // X assumed sorted and non-missing
    // w assumed non-negative and non-missing
    real scalar    i
    real colvector o
    real matrix    Q
    pointer (real matrix) scalar XX, ww
    
    if (args()<2) w = 1
    if (args()<3) P = (0, .25, .50, .75, 1)'
    if (args()<4) d = 2
    if (args()<5) fw = 0
    // handle weights
    XX = &X; ww = &w
    if (rows(w)!=1) {
        if (rows(w)!=rows(X)) _error(3200)
        // drop zero frequency observations
        if (anyof(w,0)) {
            o = select(1::rows(w), w)
            if (cols(o)==0) o = J(0,1,.) // select() may return 0x0
            XX = &(X[o,])
            ww = &(w[o])
        }
        // weights can be ignored if constant and d < 3
        if (d<3 & rows(*ww)) { // (*ww may be void)
            if (allof(*ww, (*ww)[1])) ww = &1
        }
        // weights can be ignored if w = 1 for all observations
        else if (allof(*ww,1)) ww = &1
    }
    else if (w==0) {
        XX = &(J(0,1,.))
        ww = &(J(0,1,.))
    }
    // compute weighted quantiles: requires sorted p
    if (*ww!=1) {  // not rows(*ww)!=1 !
        if (rows(P)==1) {
            o = order(P',1)
            Q = o // just to dimension q
            Q[o] = _mm_quantile_w(*XX, *ww, P[o]', d, fw, wd)
            return(Q')
        }
        Q = J(rows(P), cols(P), .)
        for (i=cols(P); i; i--) {
            o = order(P[,i],1)
            Q[o,i] = _mm_quantile_w(*XX, *ww, P[o,i], d, fw, wd)
        }
        return(Q)
    }
    // compute unweighted quantiles
    if (rows(P)==1) return(_mm_quantile_d(*XX, editmissing(P',1), d, wd)')
    Q = J(rows(P), cols(P), .)
    for (i=cols(P); i; i--) {
        Q[,i] = _mm_quantile_d(*XX, editmissing(P[,i],1), d, wd)
    }
    return(Q)
}

real colvector _mm_quantile_d(real colvector X, real colvector p, 
    real scalar d, real scalar wd)
{   // X assumed sorted and non-missing
    // p assumed nonmissing
    real scalar     n, eps
    real colvector  pn, j, j1, h

    n = rows(X)
    if ((rows(p)*n)==0) return(J(rows(p), 1, .))   // no obs or rows(p)==0
    if (n==1)           return(J(rows(p), 1, X))   // only one obs
    if (d==10) return(_mm_quantile_d_hd(X, p, wd)) // Harrell-Davis
    if (d==11) return(_mm_quantile_11(X, 1, p))    // mid-quantile by Ma et al.
    if      (d==0) pn = p * n
    else if (d==1) pn = p * n
    else if (d==2) pn = p * n
    else if (d==3) pn = p * n :- .5
                                               // pn = a + p*(n + 1 - a - b)
    else if (d==4) pn = p * n                  //      a = 0, b = 1
    else if (d==5) pn = p * n :+ .5            //      a = b = 0
    else if (d==6) pn = p * (n + 1)            //      a = b = 0 
    else if (d==7) pn = 1   :+ p * (n - 1)     //      a = b = 1
    else if (d==8) pn = 1/3 :+ p * (n + 1/3)   //      a = b = 1/3
    else if (d==9) pn = 3/8 :+ p * (n + 1/4)   //      a = b = 3/8
    else {
        display("{err}{it:def} must be an integer in [0,11]")
        _error(3300)
    }
    if (d==0) return(X[mm_clip(floor(pn):+1, 1, n)])
    if (d==1) return(X[mm_clip(ceil(pn), 1, n)])
    if (d<=3) {
        j = floor(pn)
        if (d==2) h = ((pn:>j) :+ 1) / 2
        else      h = (pn:>j) :| mod(j,2)
    }
    else {
        eps = 4 * epsilon(1) // handle rounding error as in R's quantile()
        j = floor(pn :+ eps)
        h = pn - j
        h = h :* (abs(h):>=eps)
    }
    j1 = mm_clip(j:+1, 1, n)
    _mm_clip(j, 1, n)
    return((1:-h) :* X[j] :+ h :* X[j1])
}

real colvector _mm_quantile_w(real colvector X, real colvector w, real colvector p, 
    real scalar d, real scalar fw, real scalar wd)
{   // X assumed sorted and non-missing
    // w assumed non-missing and *strictly* positive
    // p assumed sorted
    real scalar n
    real matrix W
    
    n = rows(X)
    if (n==0)       return(J(rows(p), 1, .))
    if (n==1)       return(J(rows(p), 1, X))
    if (rows(w)==0) return(J(rows(p), 1, .))
    // mid-quantile by Ma et al. (2011)
    if (d==11) return(__mm_quantile_11(X, w, p))
    // definition 3 with fw==0
    if (fw==0 & d==3) {
        if (rows(w)==1) W = (1::n) * w
        else {
            W = mm_colrunsum(w, 1, 1) // use quad precision
            if (W[n]>=.) W = J(n,1,.)
        }
        return(_mm_quantile_w_3b(X, W, p))
    }
    // other definitions
    W = _mm_ecdf2(X, w, 0, 1) // W = (uniq X, runningsum(w))
    if (fw==0 & d>3) {
        // definitions 4-10: rescale weights to effective sample size
        if (rows(w)!=1) W[,2] = W[,2] * (W[rows(W),2] / quadsum(w:^2))
    }
    if (d==0)  return(_mm_quantile_w_0(W[,1], W[,2], p))
    if (d==1)  return(_mm_quantile_w_1(W[,1], W[,2], p))
    if (d==2)  return(_mm_quantile_w_2(W[,1], W[,2], p))
    if (d==3)  return(_mm_quantile_w_3(W[,1], W[,2], p))
    if (d==4)  return(_mm_quantile_w_d(W[,1], W[,2], p, d))
    if (d==5)  return(_mm_quantile_w_d(W[,1], W[,2], p, d))
    if (d==6)  return(_mm_quantile_w_d(W[,1], W[,2], p, d))
    if (d==7)  return(_mm_quantile_w_d(W[,1], W[,2], p, d))
    if (d==8)  return(_mm_quantile_w_d(W[,1], W[,2], p, d))
    if (d==9)  return(_mm_quantile_w_d(W[,1], W[,2], p, d))
    if (d==10) return(_mm_quantile_w_hd(W[,1], W[,2], p, wd))
    display("{err}{it:def} must be an integer in [0,11]")
    _error(3300)
}

real colvector _mm_quantile_w_0(real colvector x, real colvector W, real colvector p)
{
    real scalar    i, j, pi
    real colvector P, q
    
    j = rows(W) 
    P = p * W[j]
    i = rows(p)
    q = J(i,1,.)
    for (; i; i--) {
        pi = P[i]
        for (; j>1; j--) {
            if (W[j-1]<=pi) break
        }
        q[i] = x[j]
    }
    return(q)
}

real colvector _mm_quantile_w_1(real colvector x, real colvector W, real colvector p)
{
    real scalar    i, j, pi
    real colvector P, q
    
    j = rows(W) 
    P = p * W[j]
    i = rows(p)
    q = J(i,1,.)
    for (; i; i--) {
        pi = P[i]
        for (; j>1; j--) {
            if (W[j-1]<pi) break
        }
        q[i] = x[j]
    }
    return(q)
}

real colvector _mm_quantile_w_2(real colvector x, real colvector W, real colvector p)
{
    real scalar    i, j, r, pi
    real colvector P, q
    
    j = r = rows(W) 
    P = p * W[j]
    i = rows(p)
    q = J(i,1,.)
    for (; i; i--) {
        pi = P[i]
        for (; j>1; j--) {
            if (W[j-1]<pi) break
        }
        if (W[j]==pi) {
            if (j==r) q[i] = x[j]
            else      q[i] = (x[j] + x[j+1])/2
        }
        else q[i] = x[j]
    }
    return(q)
}

real colvector _mm_quantile_w_3(real colvector x, real colvector W, real colvector p)
{
    real scalar    i, j, pi, lo, up, d0, d1
    real colvector P, q
    
    j = rows(W)
    P = p * W[j]
    i = rows(p)
    q = J(i,1,.)
    for (; i; i--) {
        pi = P[i]
        // find j such that W[j-1] < pi <= W[j]
        for (; j>1; j--) {
            if (W[j-1]<pi) break
        }
        // reached bottom
        if (j==1) { 
            q[|1\i|] = J(i,1, x[1])
            break
        }
        // if pi is larger than or equal to the next integer after W[j-1],
        // then use the upper x-value; decide between lower and upper x-value
        // only if pi is between W[j-1] and next integer
        lo = W[j-1]
        up = floor(lo) + 1
        if (pi>=up) q[i] = x[j]
        else {
            // set upper boundary to minimum between W[j] and next integer
            // after W[j-1]
            up = min( (W[j], up) )
            // obtain distances between pi and boundaries
            d0 = pi - lo
            d1 = up - pi
            // case 1: pi closer to upper boundary
            if (d1<d0) q[i] = x[j]
            // case 2: pi closer to lower boundary
            else if (d0<d1) q[i] = x[j-1]
            // case 3: equal distance
            else {
                // use x-value from lower boundary if lower boundary is integer 
                // and even; else use x-value from upper boundary; this is an
                // arbitrary decision
                if (!mod(lo,2)) q[i] = x[j-1]
                else            q[i] = x[j]
            }
        }
    }
    return(q)
}

real colvector _mm_quantile_w_3b(real colvector x, real colvector W, real colvector p)
{
    real scalar    i, j, pi, d0, d1
    real colvector P, q
    
    j = rows(W)
    P = p * W[j]
    i = rows(p)
    q = J(i,1,.)
    for (; i; i--) {
        pi = P[i]
        for (; j>1; j--) {
            if (W[j-1]<pi) break
        }
        // reached bottom
        if (j==1) { 
            q[|1\i|] = J(i,1, x[1])
            break
        }
        // obtain distances
        d0 = pi - W[j-1]
        d1 = W[j] - pi
        // case 1: pi closer to upper boundary
        if (d1<d0) q[i] = x[j]
        // case 2: pi closer to lower boundary
        else if (d0<d1) q[i] = x[j-1]
        // case 3: equal distance
        else {
            // use x-value from upper boundary if index of upper boundary is
            // even; else use x-value from lower boundary
            if (!mod(j,2)) q[i] = x[j]
            else           q[i] = x[j-1]
        }
    }
    return(q)
}

real colvector _mm_quantile_w_d(real colvector x, real colvector W, 
    real colvector p, real scalar d) // modifies W
{
    real scalar    N, n, i, r, ll, ul, a, b, hmax
    real colvector h, q, ww
    
    N = rows(x)
    n = W[N]
    if      (d==4) h = p * n
    else if (d==5) h = p * n :+ .5
    else if (d==6) h = p * (n + 1)
    else if (d==7) h = 1   :+ p * (n - 1)
    else if (d==8) h = 1/3 :+ p * (n + 1/3)
    else if (d==9) h = 3/8 :+ p * (n + 1/4)
    r = rows(h)
    if (r) {
        if ((h[1]-1)<=W[1]) {
            // duplicate min if lower limit of first quantile is below data
            x = x[1] \ x 
            W = (h[1]-2) \ W
            N = N + 1
        }
        if (h[r]>W[N]) {
            // duplicate max if upper limit of last quantile is above data
            // (using max(h) instead of h[r] because p may include missing)
            hmax = max(h)
            if (hmax<. & hmax>W[N]) {
                x = x \ x[N]
                W = W \ hmax
                N = N + 1
            }
        }
    }
    q = J(r,1,.)
    a = 2
    for (i=1; i<=r; i++) {
        ul = h[i]; ll = ul - 1
        // exit once lower limit is above data (this also takes care of p>=.)
        if (ll>W[N]) {
            q[|i\.|] = J(r-i+1, 1, x[N])
            break
        }
        // lower index
        while (1) {
            if (W[a]>=ll) break
            if (a>=N) break
            a++
        }
        // upper upper index
        b = a
        if (b<N) {
            while (1) {
                if (W[b]>=ul) break
                if (b>=N) break
                b++
            }
        }
        // compute quantile
        ww = W[|a-1 \ b|] :- ll
        ww = mm_diff(rowmax((J(rows(ww),1,0), rowmin((ww, J(rows(ww),1,1))))))
        q[i] = quadsum(ww :* x[|a \ b|])
    }
    return(q)
}

real colvector _mm_quantile_d_hd(real colvector X, real colvector p,
    real scalar wd)
{
    real scalar    n, i
    real colvector q, W
    
    i = rows(p)
    q = J(i,1,.)
    n = rows(X)
    W = 0 \ (1::n)/n
    for (; i; i--) q[i] = __mm_quantile_hd(X, W, p[i], n, wd)
    return(q)
}

real colvector _mm_quantile_w_hd(real colvector X, real colvector W,
    real colvector p, real scalar wd)  // modifies W
{
    real scalar    n, i
    real colvector q
    
    i = rows(p)
    q = J(i,1,.)
    n = W[rows(W)] 
    W = 0 \ W/n
    for (; i; i--) q[i] = __mm_quantile_hd(X, W, max((0,min((1,p[i])))), n, wd)
    return(q)
}

real scalar __mm_quantile_hd(real colvector X, real colvector W, real scalar p,
    real scalar n, real scalar wd0)
{   // note: domain of a and b in ibeta() is 1e-10 to 1e+17; see help ibeta()
    real scalar    a, b, wd, l, r, il, ir
    real colvector w
    
    a = (n + 1) * p
    b = (n + 1) * (1 - p)
    if      (a<1e-10) return(X[1])
    else if (b<1e-10) return(X[rows(X)])
    if (a>1e+17 | b>1e+17) {
        display("{err}sample size too large; cannot evaluate {bf:ibeta()}")
        _error(3498)
    }
    // trimmed estimator (Akinshin 2021)
    if (wd0<1) {
        wd = wd0<=0 ? 1 / sqrt(n) : wd0
        if      (a<=1 & b<=1) {; l=.   ; r=. ; } // can only happen if n=1
        else if (a<=1 & b>1 ) {; l=0   ; r=wd; } // left border
        else if (a>1  & b<=1) {; l=1-wd; r=1 ; } // right border
        else {; l = __mm_quantile_hd_l(a, b, wd); r = l+wd; }
        if (l<.) { // (else use untrimmed estimator)
            // find lower index
            if (l<=0) il = 1
            else      mm_hunt(W, l, il = floor(l*rows(X)))
            // find upper index
            if (r>=1) ir = rows(W)
            else      {; mm_hunt(W, r, ir = il); ir++; }
            // compute weights
            w = W[|il\ir|]; w[1] = l; w[rows(w)] = r
            w = (ibeta(a, b, w) :- ibeta(a,b,l)) / mm_diff(ibeta(a, b, l\r))
            return(quadsum(mm_diff(w) :* X[|il\ir-1|]))
        }
    }
    // untrimmed estimator
    w = mm_diff(ibeta(a, b, W))
    return(quadsum(w :* X)) // (faster than quadcross(w, X))
}

real scalar __mm_quantile_hd_l(real scalar a, real scalar b, real scalar wd)
{
    real scalar m, lo, up, l, rc
    
    m = (a - 1) / (a + b - 2)
    lo = max((0, m-wd))
    up = min((m, 1-wd))
    rc = mm_root(l=., &___mm_quantile_hd_l(), lo, up, 0, 1000, a, b, wd)
    if (rc) _error(3498, "could not locate highest density interval")
    return(l)
}

real scalar ___mm_quantile_hd_l(real scalar x, real scalar a, real scalar b,
    real scalar wd) return(betaden(a, b, x) - betaden(a, b, x+wd))

real colvector _mm_quantile_11(real colvector X, real colvector w,
    real colvector p)
{   // mid-quantile by Ma et al. (2011)
    // X assumed sorted and non-missing
    // w assumed non-missing and *strictly* positive
    // p not assumed sorted
    real colvector o
    real colvector q
    
    o = order(p,1)
    q = o   // just to dimension q
    q[o] = __mm_quantile_11(X, w, p[o])
    return(q)
}

real colvector __mm_quantile_11(real colvector X, real colvector w,
    real colvector p)
{   // mid-quantile by Ma et al. (2011)
    // X assumed sorted and non-missing
    // w assumed non-missing and *strictly* positive
    // p assumed sorted
    real scalar    i, pi, j, xj, pj
    real colvector q
    real matrix    m
    
    i = rows(p)
    q = J(i,1,.)
    m = _mm_ecdf2(X, w, 1) // mid cdf at unique values: (x, mid cdf)
    j = rows(m)
    // handle p > max(cdf)
    xj = m[j,1]; pj = m[j,2]
    for (; i; i--) {
        if (pj>=p[i]) break
        q[i] = xj
    }
    // handle rest
    for (; i; i--) {
        pi = p[i]
        while (j) {
            pj = m[j,2]
            if (pj<pi) break
            j--
        }
        if (!j) { // reached bottom (p <= min(cdf))
            q[|1\i|] = J(i, 1, m[1,1])
            break
        }
        // note: j must be lower than rows(m) at this point
        xj = m[j,1]
        q[i] = xj + (pi-pj)/(m[j+1,2]-pj) * (m[j+1,1] - xj)
    }
    return(q)
}

end
