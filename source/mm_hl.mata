*! version 1.0.0  21oct2020  Ben Jann
version 11
mata:

// robust statistics based on pairwise comparisons (code adapted from
// robstat.ado, version 1.0.3)
// - HL: location
// - Qn: scale
// - mc (medcouple): skewness

// HL estimator

real scalar mm_hl(real colvector X, | real colvector w, real scalar fw, 
    real scalar naive)
{
    real colvector p
    
    if (args()<2) w = 1
    if (args()<3) fw = 0
    if (args()<4) naive = 0
    if (hasmissing(w)) _error(3351)
    if (any(w:<0)) _error(3498, "negative weights not allowed")
    if (rows(X)==0) return(.)
    if (hasmissing(X)) _error(3351)
    if (naive) {
        if (rows(w)==1) return(_mm_hl_naive(X))
        if (fw) return(_mm_hl_naive_fw(X, w))
        return(_mm_hl_naive_w(X, w))
    }
    if (mm_issorted(X)) {
        if (rows(w)==1) return(_mm_hl(X))
        if (fw) return(_mm_hl_fw(X, w))
        return(_mm_hl_w(X, w))
    }
    if (rows(w)==1) return(_mm_hl(sort(X,1)))
    p = order((X,w), (1,2))
    if (fw) return(_mm_hl_fw(X[p], w[p]))
    return(_mm_hl_w(X[p], w[p]))
}

real scalar _mm_hl(real colvector x) // no weights
{
    // the trick of this algorithm is to consider only elements that are on the 
    // right of the main diagonal
    real scalar     i, j, m, k, n, nl, nr, nL, nR, trial
    real colvector  xx, /*ww,*/ l, r, L, R
    
    n = rows(x)                    // dimension of search matrix
    if (n==1) return(x)            // returning observed value if n=1
    xx      = /*ww =*/ J(n, 1, .)  // temp vector for matrix elements
    l = L   = (1::n):+1            // indices of left boundary (old and new)
    r = R   = J(n, 1, n)           // indices of right boundary (old and new)
    nl = nl = n + comb(n, 2)       // number of cells below left boundary
    nr = nR = n * n                // number of cells within right boundary
    k       = nl + comb(n, 2)/2    // target quantile
    while ((nr-nl)>n) {
        // get trial value
        m = 0
        for (i=1; i<n; i++) { // last row cannot contain candidates
            if (l[i]<=r[i]) {
                // high median within row
                xx[++m] = __mm_hl_el(x, i, l[i]+trunc((r[i]-l[i]+1)/2))
                /*m++
                ww[m] = r[i] - l[i] + 1
                xx[m] = __mm_hl_el(x, i, l[i]+trunc(ww[m]/2))*/
            }
        }
        trial = _mm_hl_qhi(xx[|1 \ m|], .5)
        /*trial = _mm_hl_qhi_w(xx[|1 \ m|], ww[|1 \ m|], .5)*/
        /*the unweighted quantile is faster; results are the same*/
        // move right border
        j = n-1
        for (i=(n-1); i>=1; i--) {
            if (i==j) {
                if (__mm_hl_el(x, i, j)>=trial) {
                    R[i] = j
                    j = i-1
                    continue
                }
            }
            if (j<n) {
                while (__mm_hl_el(x, i, j+1)<trial) {
                    j++
                    if (j==n) break
                }
            }
            R[i] = j
        }
        nR = sum(R)
        if (nR>k) {
            swap(r, R)
            nr = nR
            continue
        }
        // move left border
        j = n + 1
        for (i=1; i<=n; i++) {
            if (j>(i+1)) {
                while (__mm_hl_el(x, i, j-1)>trial) {
                    j--
                    if (j==(i+1)) break
                }
            }
            if (j<(i+1)) j = (i+1)
            L[i] = j
        }
        nL = sum(L) - n
        if (nL<k) {
            swap(l, L)
            nl = nL
            continue
        }
        // trial = low quantile = high quantile
        if (ceil(k)!=k | (nR<k & nL>k)) return(trial)
        // trial = low quantile
        if (nL==k) {
            m = 0
            for (i=1; i<=n; i++) {
                if (L[i]>n) continue
                xx[++m] = __mm_hl_el(x, i, L[i])
            }
            return((trial+min(xx[|1 \ m|]))/2)
        }
        // trial = high quantile
        m = 0
        for (i=1; i<=n; i++) {
            if (R[i]<=i) continue
            xx[++m] = __mm_hl_el(x, i, R[i])
        }
        return((trial+max(xx[|1 \ m|]))/2)
    }
    // get target value from remaining candidates
    m = 0
    for (i=1; i<n; i++) { // last row cannot contain candidates
        if (l[i]<=r[i]) {
            for (j=l[i]; j<=r[i]; j++) {
                m++
                xx[m] = __mm_hl_el(x, i, j)
            }
        }
    }
    return(_mm_hl_q(xx[|1 \ m|], k, nl))
}

real scalar __mm_hl_el(real colvector y, real scalar i, real scalar j)
{
    return((y[i] + y[j])/2)
}

real scalar _mm_hl_w(real colvector x, real colvector w)
{
    real scalar     i, j, m, k, n, nl, nr, trial, Wl, WR, WL, W0, W1
    real colvector  xx, ww, l, r, L, R, ccw
    
    n       = rows(x)              // dimension of search matrix
    if (n==1) return(x)            // returning observed value if n=1
    xx = ww = J(n, 1, .)           // temp vector for matrix elements
    l = L   = (1::n):+1            // indices of left boundary (old and new)
    r = R   = J(n, 1, n)           // indices of right boundary (old and new)
    nl      = comb(n, 2) + n       // number of cells below left boundary
    nr      = n * n                // number of cells within right boundary
    ccw     = quadrunningsum(w)    // cumulative column weights
    W1      = quadsum(w[|2 \ .|] :* ccw[|1 \ n-1|]) // sum of weights in target triangle
    W0      = W1 + quadsum(w:*w)   // sum of weights in rest of search matrix
    Wl = WL = W0                   // sum of weights below left boundary
    WR      = W0 + W1              // sum of weights within right boundary
    k       = W0 + W1/2            // target quantile
    while ((nr-nl)>n) {
        // get trial value
        m = 0
        for (i=1; i<n; i++) { // last row cannot contain candidates
            if (l[i]<=r[i]) {
                // high median within row
                xx[++m] = __mm_hl_el(x, i, l[i]+trunc((r[i]-l[i]+1)/2))
            }
        }
        trial = _mm_hl_qhi(xx[|1 \ m|], .5)
        // move right border
        j = n-1
        for (i=(n-1); i>=1; i--) {
            if (i==j) {
                if (__mm_hl_el(x, i, j)>=trial) {
                    R[i] = j
                    j = i-1
                    continue
                }
            }
            if (j<n) {
                while (__mm_hl_el(x, i, j+1)<trial) {
                    j++
                    if (j==n) break
                }
            }
            R[i] = j
        }
        WR = quadsum(w:*ccw[R])
        if (WR>k) {
            swap(r, R)
            nr = sum(R)
            continue
        }
        // move left border
        j = n + 1
        for (i=1; i<=n; i++) {
            if (j>(i+1)) {
                while (__mm_hl_el(x, i, j-1)>trial) {
                    j--
                    if (j==(i+1)) break
                }
            }
            if (j<(i+1)) j = (i+1)
            L[i] = j
        }
        WL = quadsum(w:*ccw[L:-1])
        if (WL<k) {
            swap(l, L)
            Wl = WL
            nl = sum(L) - n
            continue
        }
        // trial = low quantile = high quantile
        if (WR==WL | (WR<k & WL>k)) return(trial)
        // trial = low quantile
        if (WL==k) {
            m = 0
            for (i=1; i<=n; i++) {
                if (L[i]>n) continue
                xx[++m] = __mm_hl_el(x, i, L[i])
            }
            return((trial+min(xx[|1 \ m|]))/2)
        }
        // trial = high quantile
        m = 0
        for (i=1; i<=n; i++) {
            if (R[i]<=i) continue
            xx[++m] = __mm_hl_el(x, i, R[i])
        }
        return((trial+max(xx[|1 \ m|]))/2)
    }
    // get target value from remaining candidates
    m = 0
    for (i=1; i<n; i++) { // last row cannot contain candidates
        if (l[i]<=r[i]) {
            for (j=l[i]; j<=r[i]; j++) {
                m++
                xx[m] = __mm_hl_el(x, i, j)
                ww[m] = w[i] * w[j]
            }
        }
    }
    return(_mm_hl_q_w(xx[|1 \ m|], ww[|1 \ m|], k, Wl))
}

real scalar _mm_hl_fw(real colvector x, real colvector w)
{   // the algorithm "duplicates" the diagonal so that the relevant pairs in 
    // case of w>1 can be taken into account
    real scalar     i, j, m, k, n, nl, nr, trial, Wl, WR, WL, W0, W1
    real colvector  xx, ww, l, r, L, R, ccw, wcorr, idx
    
    if (any(trunc(w):!=w)) _error(3498, "non-integer frequency not allowed")
    n       = rows(x)          // dimension of search matrix
    if (n==1) return(x)        // returning observed value if n=1
    xx = ww = J(n, 1, .)       // temp vector for matrix elements
    ccw     = runningsum(w)    // cumulative column weights
    wcorr   = mm_cond(w:<=1, 0, comb(w, 2)) // correction of weights
    idx     = 1::n             // diagonal indices
    l = L   = idx :+ 1 :+ (wcorr:==0) // indices of left boundary (old and new)
    r = R   = J(n, 1, n+1)     // indices of right boundary (old and new)
    nl      = comb(n, 2) + n + sum(wcorr:==0) // n. of cells below left boundary
    nr      = n * (n+1)        // number of cells within right boundary
    W0 = W1 = sum(w[|2 \ .|] :* ccw[|1 \ n-1|])
    W1      = W1 + sum(wcorr)  // sum of weights in target triangle
    W0      = W0 + sum(w) + sum(wcorr) // sum of weights in rest of search matrix
    Wl = WL = W0               // sum of weights below left boundary
    WR      = W0 + W1          // sum of weights within right boundary
    k       = W0 + W1/2        // target quantile
    while ((nr-nl)>n) {
        // get trial value
        m = 0
        for (i=1; i<=n; i++) {
            if (l[i]<=r[i]) {
                // high median within row
                xx[++m] = __mm_hl_el(x, i, 
                    (l[i]-1) + trunc(((r[i]-1)-(l[i]-1)+1)/2))
            }
        }
        trial = _mm_hl_qhi(xx[|1 \ m|], .5)
        // move right border
        j = n
        for (i=n; i>=1; i--) {
            if (j==i) {
                if (__mm_hl_el(x, i, j)>=trial) {
                    R[i] = j + (wcorr[i]==0)
                    j = i-1
                    continue
                }
            }
            if (j<=n) {
                while (__mm_hl_el(x, i, (j+1)-((j+1)>i))<trial) {
                    j++
                    if (j>n) break
                }
            }
            R[i] = j
        }
        WR = sum(w:*ccw[R:-(R:>idx)]) - sum(wcorr:*(R:==idx))
        if (WR>k) {
            swap(r, R)
            nr = sum(R)
            continue
        }
        // move left border
        j = n + 2
        for (i=1; i<=n; i++) {
            if (j>(i+1)) {
                while (__mm_hl_el(x, i, (j-1)-((j-1)>i))>trial) {
                    j--
                    if (j==(i+1)) break
                }
            }
            if (j<(i+1)) j = (i+1)
            if (j==(i+1)) {
                if (wcorr[i]==0) j++
            }
            L[i] = j
        }
        WL = sum(w:*ccw[L:-1:-((L:-1):>idx)]) - sum(wcorr:*((L:-1):==idx))
        if (WL<k) {
            swap(l, L)
            Wl = WL
            nl = sum(L) - n
            continue
        }
        // trial = low quantile = high quantile
        if (ceil(k)!=k | (WR<k & WL>k)) return(trial)
        // trial = low quantile
        if (WL==k) {
            m = 0
            for (i=1; i<=n; i++) {
                if ((L[i]-(L[i]>i))>n) continue
                xx[++m] = __mm_hl_el(x, i, L[i]-(L[i]>i))
            }
            return((trial+min(xx[|1 \ m|]))/2)
        }
        // trial = high quantile
        m = 0
        for (i=1; i<=n; i++) {
            if (R[i]<=i) continue
            if (R[i]==i+1) {
                if (wcorr[i]==0) continue
            }
            xx[++m] = __mm_hl_el(x, i, R[i]-1)
        }
        return((trial+max(xx[|1 \ m|]))/2)
    }
    // get target value from remaining candidates
    m = 0
    for (i=1; i<=n; i++) {
        if (l[i]<=r[i]) {
            for (j=l[i]; j<=r[i]; j++) {
                m++
                xx[m] = __mm_hl_el(x, i, j-1)
                ww[m] = (i==(j-1) ? comb(w[i], 2) : w[i]*w[j-1])
            }
        }
    }
    return(_mm_hl_q_w(xx[|1 \ m|], ww[|1 \ m|], k, Wl))
}

real scalar _mm_hl_naive(real colvector x) // no weights
{
    real scalar    i, j, m, n
    real colvector xx
    
    n = rows(x)
    if (n==1) return(x) // HL undefined if n=1; returning observed value
    m = 0
    xx = J(comb(n,2), 1, .)
    for (i=1; i<n; i++) {
        for (j=(i+1); j<=n; j++) {
            xx[++m] = __mm_hl_el(x, i, j)
        }
    }
    return(_mm_hl_q(xx, .5))
}

real scalar _mm_hl_naive_w(real colvector x, real colvector w)
{
    real scalar    i, j, m, n
    real colvector xx, ww

    n = rows(x)
    if (n==1) return(x) // HL undefined if n=1; returning observed value
    m = 0
    xx = J(comb(n,2), 1, .)
    ww = J(rows(xx), 1, .)
    for (i=1; i<n; i++) {
        for (j=(i+1); j<=n; j++) {
            m++
            xx[m] = __mm_hl_el(x, i, j)
            ww[m] = w[i]*w[j]
        }
    }
    return(_mm_hl_q_w(xx, ww, .5))
}

real scalar _mm_hl_naive_fw(real colvector x, real colvector w)
{
    real scalar    i, j, m, n
    real colvector xx, ww

    if (any(trunc(w):!=w)) _error(3498, "non-integer frequency not allowed")
    n = rows(x)
    if (n==1) return(x) // HL undefined if n=1; returning observed value
    m = 0
    xx = J(comb(n,2)+sum(w:>1), 1, .)
    ww = J(rows(xx), 1, .)
    for (i=1; i<=n; i++) {
        for (j=i; j<=n; j++) {
            if (i==j) {
                if (w[i]==1) continue
            }
            m++
            xx[m] = __mm_hl_el(x, i, j)
            ww[m] = (i==j ? comb(w[i], 2) : w[i]*w[j])
        }
    }
    return(_mm_hl_q_w(xx, ww, .5))
}

// Qn estimator

real scalar mm_qn(real colvector X, | real colvector w, real scalar fw, 
    real scalar naive)
{
    real colvector p
    real scalar    c
    
    if (args()<2) w = 1
    if (args()<3) fw = 0
    if (args()<4) naive = 0
    if (hasmissing(w)) _error(3351)
    if (any(w:<0)) _error(3498, "negative weights not allowed")
    if (rows(X)==0) return(.)
    if (hasmissing(X)) _error(3351)
    c = 1 / (sqrt(2) * invnormal(5/8))
    if (naive) {
        if (rows(w)==1) return(_mm_qn_naive(X) * c)
        if (fw) return(_mm_qn_naive_fw(X, w) * c)
        return(_mm_qn_naive_w(X, w) * c)
    }
    if (mm_issorted(X)) {
        if (rows(w)==1) return(_mm_qn(X) * c)
        if (fw) return(_mm_qn_fw(X, w) * c)
        return(_mm_qn_w(X, w) * c)
    }
    if (rows(w)==1) return(_mm_qn(sort(X,1)) * c)
    p = order((X,w), (1,2))
    if (fw) return(_mm_qn_fw(X[p], w[p]) * c)
    return(_mm_qn_w(X[p], w[p]) * c)
}

real scalar _mm_qn(real colvector x) // no weights
{
    real scalar     i, j, m, k, n, nl, nr, nL, nR, trial
    real colvector  xx, /*ww,*/ l, r, L, R

    n = rows(x)                    // dimension of search matrix
    if (n==1) return(0)            // returning zero if n=1
    xx      = /*ww =*/ J(n, 1, .)  // temp vector for matrix elements
    l = L   = (n::1):+1            // indices of left boundary (old and new)
    r = R   = J(n, 1, n)           // indices of right boundary (old and new)
    nl = nl = comb(n, 2) + n       // number of cells below left boundary
    nr = nR = n * n                // number of cells within right boundary
    k       = nl + comb(n, 2)/4    // target quantile
    /*k = nl + comb(trunc(n/2) + 1, 2)*/
    while ((nr-nl)>n) {
        // get trial value
        m = 0
        for (i=2; i<=n; i++) { // first row cannot contain candidates
            if (l[i]<=r[i]) {
                // high median within row
                xx[++m] = __mm_qn_el(x, i, l[i]+trunc((r[i]-l[i]+1)/2), n)
                /*m++
                ww[m] = r[i] - l[i] + 1
                xx[m] = __mm_qn_el(x, i, l[i]+trunc(www[m]/2), n)*/
            }
        }
        trial = _mm_hl_qhi(xx[|1 \ m|], .5)
        /*trial = _mm_hl_qhi_w(xx[|1 \ m|], ww[|1 \ m|], .5)*/
        /*the unweighted quantile is faster; results are the same*/
        //move right border
        j = 0
        for (i=n; i>=1; i--) {
            if (j<n) {
                while (__mm_qn_el(x, i, j+1, n)<trial) {
                    j++
                    if (j==n) break
                }
            }
            R[i] = j
        }
        nR = sum(R)
        if (nR>k) {
            swap(r, R)
            nr = nR
            continue
        }
        // move left border
        j = n + 1
        for (i=1; i<=n; i++) {
            while (__mm_qn_el(x, i, j-1, n)>trial) {
                j--
            }
            L[i] = j
        }
        nL = sum(L) - n
        if (nL<k) {
            swap(l, L)
            nl = nL
            continue
        }
        // trial = low quantile = high quantile
        if (ceil(k)!=k | (nR<k & nL>k)) return(trial)
        // trial = low quantile
        if (nL==k) {
            m = 0
            for (i=1; i<=n; i++) {
                if (L[i]>n) continue
                xx[++m] = __mm_qn_el(x, i, L[i], n)
            }
            return((trial+min(xx[|1 \ m|]))/2)
        }
        // trial = high quantile
        for (i=1; i<=n; i++) {
            xx[i] = __mm_qn_el(x, i, R[i], n)
        }
        return((trial+max(xx))/2)
    }
    // get target value from remaining candidates
    m = 0
    for (i=2; i<=n; i++) { // first row cannot contain candidates
        if (l[i]<=r[i]) {
            for (j=l[i]; j<=r[i]; j++) {
                m++
                xx[m] = __mm_qn_el(x, i, j, n)
            }
        }
    }
    return(_mm_hl_q(xx[|1 \ m|], k, nl))
}

real scalar __mm_qn_el(real colvector y, real scalar i, real scalar j, 
    real scalar n)
{
    return(y[i] - y[n-j+1])
}

real scalar _mm_qn_w(real colvector x, real colvector w)
{
    real scalar     i, j, m, k, n, nl, nr, trial, Wl, WR, WL, W0, W1
    real colvector  xx, ww, l, r, L, R, p, ccw

    n       = rows(x)              // dimension of search matrix
    if (n==1) return(0)            // returning zero if n=1
    xx = ww = J(n, 1, .)           // temp vector for matrix elements
    l = L   = (n::1):+1            // indices of left boundary (old and new)
    r = R   = J(n, 1, n)           // indices of right boundary (old and new)
    nl      = comb(n, 2) + n       // number of cells below left boundary
    nr      = n * n                // number of cells within right boundary
    ccw     = quadrunningsum(w[n::1]) // cumulative column weights
    W0      = quadsum(w:*ccw[n::1]) // sum weights in rest of search matrix
    W1      = quadsum(w:*ccw[rows(ccw)]) - W0 // sum of weights target triangle
    Wl = WL = W0                   // sum of weights below left boundary
    WR      = W0 + W1              // sum of weights within right boundary
    k       = W0 + W1/4            // target sum (high 25% quantile)
    while ((nr-nl)>n) {
        // get trial value
        m = 0
        for (i=2; i<=n; i++) { // first row cannot contain candidates
            if (l[i]<=r[i]) {
                // high median within row
                xx[++m] = __mm_qn_el(x, i, l[i]+trunc((r[i]-l[i]+1)/2), n)
            }
        }
        trial = _mm_hl_qhi(xx[|1 \ m|], .5)
        //move right border
        j = 0
        for (i=n; i>=1; i--) {
            if (j<n) {
                while (__mm_qn_el(x, i, j+1, n)<trial) {
                    j++
                    if (j==n) break
                }
            }
            R[i] = j
        }
        p = (R:>0)
        WR = quadsum(select(w, p) :* ccw[select(R, p)])
        if (WR>k) {
            swap(r, R)
            nr = sum(R)
            continue
        }
        // move left border
        j = n + 1
        for (i=1; i<=n; i++) {
            while (__mm_qn_el(x, i, j-1, n)>trial) {
                j--
            }
            L[i] = j
        }
        WL = quadsum(w :* ccw[L:-1])
        if (WL<k) {
            swap(l, L)
            Wl = WL
            nl = sum(L) - n
            continue
        }
        // trial = low quantile = high quantile
        if (WR==WL | (WR<k & WL>k)) return(trial)
        // trial = low quantile
        if (WL==k) {
            m = 0
            for (i=1; i<=n; i++) {
                if (L[i]>n) continue
                xx[++m] =  __mm_qn_el(x, i, L[i], n)
            }
            return((trial+min(xx[|1 \ m|]))/2)
        }
        // trial = high quantile
        for (i=1; i<=n; i++) {
            xx[i] = __mm_qn_el(x, i, R[i], n)
        }
        return((trial+max(xx))/2)
    }
    // get target value from remaining candidates
    m = 0
    for (i=2; i<=n; i++) { // first row cannot contain candidates
        if (l[i]<=r[i]) {
            for (j=l[i]; j<=r[i]; j++) {
                m++
                xx[m] = __mm_qn_el(x, i, j, n)
                ww[m] = w[i] * w[n-j+1]
            }
        }
    }
    return(_mm_hl_q_w(xx[|1 \ m|], ww[|1 \ m|], k, Wl))
}

real scalar _mm_qn_fw(real colvector x, real colvector w)
{
    real scalar     i, j, m, k, n, nl, nr, trial, Wl, WR, WL, W0, W1
    real colvector  xx, ww, l, r, L, R, p, ccw, wcorr, idx

    if (any(trunc(w):!=w)) _error(3498, "non-integer frequency not allowed")
    n       = rows(x)              // dimension of search matrix
    if (n==1) return(0)            // returning zero if n=1
    xx = ww = J(n, 1, .)           // temp vector for matrix elements
    idx     = n::1                 // (minor) diagonal indices
    l = L   = idx:+1               // indices of left boundary (old and new)
    r = R   = J(n, 1, n+1)         // indices of right boundary (old and new)
    nl      = comb(n, 2) + n       // number of cells below left boundary
    nr      = n * (n+1)            // number of cells within right boundary
    ccw     = runningsum(w[idx])   // cumulative column weights
    wcorr   = mm_cond(w:<=1, 0, comb(w, 2))[idx] // correction of weights
    W0      = sum(w:*ccw[idx]) - sum(wcorr) // sum of weights in rest of search matrix
    W1      = sum(w:*ccw[rows(ccw)]) - W0 // sum of weights in target triangle
    Wl = WL = W0                   // sum of weights below left boundary
    WR      = W0 + W1              // sum of weights within right boundary
    k       = W0 + W1/4            // target quantile
    while ((nr-nl)>n) {
        // get trial value
        m = 0
        for (i=1; i<=n; i++) {
            if (l[i]<=r[i]) {
                // high median within row
                xx[++m] = __mm_qn_el(x, i, 
                    (l[i]-1)+trunc(((r[i]-1)-(l[i]-1)+1)/2), n)
            }
        }
        trial = _mm_hl_qhi(xx[|1 \ m|], .5)
        //move right border
        j = 0
        for (i=n; i>=1; i--) {
            if (j<=n) {
                while (__mm_qn_el(x, i, (j+1)-((j+1)>(n-i+1)), n)<trial) {
                    j++
                    if (j>n) break
                }
            }
            R[i] = j
        }
        p = (R:>0)
        WR = sum(select(w, p) :* ccw[select(R:-(R:>idx), p)]) - sum(wcorr:*(R:==idx))
        if (WR>k) {
            swap(r, R)
            nr = sum(R)
            continue
        }
        // move left border
        j = n + 2
        for (i=1; i<=n; i++) {
            while (__mm_qn_el(x, i, (j-1)-((j-1)>(n-i+1)), n)>trial) {
                j--
            }
            L[i] = j
        }
        WL = sum(w:*ccw[L:-1:-((L:-1):>idx)]) - sum(wcorr:*((L:-1):==idx))
        if (WL<k) {
            swap(l, L)
            Wl = WL
            nl = sum(L) - n
            continue
        }
        // trial = low quantile = high quantile
        if (ceil(k)!=k | (WR<k & WL>k)) return(trial)
        // trial = low quantile
        if (WL==k) {
            m = 0
            for (i=1; i<=n; i++) {
                if ((L[i]-(L[i]>(n-i+1)))>n) continue
                xx[++m] = __mm_qn_el(x, i, L[i]-(L[i]>(n-i+1)), n)
            }
            return((trial+min(xx[|1 \ m|]))/2)
        }
        // trial = high quantile
        for (i=1; i<=n; i++) {
            xx[i] = __mm_qn_el(x, i, R[i]-(R[i]>(n-i+1)), n)
        }
        return((trial+max(xx))/2)
    }
    // get target value from remaining candidates
    m = 0
    for (i=1; i<=n; i++) {
        if (l[i]<=r[i]) {
            for (j=l[i]; j<=r[i]; j++) {
                m++
                xx[m] = __mm_qn_el(x, i, j-1, n)
                ww[m] = ((n-i+1)==(j-1) ? comb(w[i], 2) : w[i]*w[n-(j-1)+1])
            }
        }
    }
    return(_mm_hl_q_w(xx[|1 \ m|], ww[|1 \ m|], k, Wl))
}

real scalar _mm_qn_naive(real colvector x) // no weights
{
    real scalar    i, j, m, n
    real colvector xx
    
    n = rows(x)
    if (n==1) return(0) // returning zero if n=1
    m = 0
    xx = J(comb(n,2), 1, .)
    for (i=1; i<n; i++) {
        for (j=(i+1); j<=n; j++) {
            xx[++m] = abs(x[i] - x[j])
        }
    }
    return(_mm_hl_q(xx, 0.25))
}

real scalar _mm_qn_naive_w(real colvector x, real colvector w)
{
    real scalar    i, j, m, n
    real colvector xx, ww
    
    n = rows(x)
    if (n==1) return(0) // Qn undefined if n=1; returning zero
    m = 0
    xx = J(comb(n,2), 1, .)
    ww = J(rows(xx), 1, .)
    for (i=1; i<n; i++) {
        for (j=(i+1); j<=n; j++) {
            m++
            xx[m] = abs(x[i] - x[j])
            ww[m] = w[i]*w[j]
        }
    }
    return(_mm_hl_q_w(xx, ww, 0.25))
}

real scalar _mm_qn_naive_fw(real colvector x, real colvector w)
{
    real scalar    i, j, m, n
    real colvector xx, ww

    if (any(trunc(w):!=w)) _error(3498, "non-integer frequency not allowed")
    n = rows(x)
    if (n==1) return(0) // returning zero if n=1
    m = 0
    xx = J(comb(n,2)+sum(w:>1), 1, .)
    ww = J(rows(xx), 1, .)
    for (i=1; i<=n; i++) {
        for (j=i; j<=n; j++) {
            if (i==j) {
                if (w[i]==1) continue
            }
            m++
            xx[m] = abs(x[i] - x[j])
            ww[m] = (i==j ? comb(w[i], 2) : w[i]*w[j])
        }
    }
    return(_mm_hl_q_w(xx, ww, 0.25))
}

// medcouple

real scalar mm_mc(real colvector X, | real colvector w, real scalar fw, 
    real scalar naive)
{
    real colvector p, r
    
    if (args()<2) w = 1
    if (args()<3) fw = 0
    if (args()<4) naive = 0
    if (hasmissing(w)) _error(3351)
    if (any(w:<0)) _error(3498, "negative weights not allowed")
    if (rows(X)==0) return(.)
    if (hasmissing(X)) _error(3351)
    if (mm_issorted(X)) {
        if (naive) {
            if (rows(w)==1) return(_mm_mc_naive(X:-_mm_median(X)))
            if (fw) return(_mm_mc_naive_fw(X:-_mm_median(X,w), w))
            return(_mm_mc_naive_w(X:-_mm_median(X,w), w))
        }
        r = rows(X)::1
        if (rows(w)==1) return(_mm_mc(X[r]:-_mm_median(X)))
        if (fw) return(_mm_mc_fw(X[r]:-_mm_median(X,w), w[r]))
        return(_mm_mc_w(X[r]:-_mm_median(X,w), w[r]))
    }
    if (naive) {
        if (rows(w)==1) return(_mm_mc_naive(X:-mm_median(X)))
        if (fw) return(_mm_mc_naive_fw(X:-mm_median(X,w), w))
        return(_mm_mc_naive_w(X:-mm_median(X,w), w))
    }
    if (rows(w)==1) {
        p = order(X, 1)
        r = p[rows(X)::1]
        return(_mm_mc(X[r]:-_mm_median(X[p])))
    }
    p = order((X,w), (1,2))
    r = p[rows(X)::1]
    if (fw) return(_mm_mc_fw(X[r]:-_mm_median(X[p], w[p]), w[r]))
    return(_mm_mc_w(X[r]:-_mm_median(X[p], w[p]), w[r]))
}

real scalar _mm_mc(real colvector x) // no weights; assumes med(x)=0
{
    real scalar     i, j, m, k, n, q, nl, nr, nL, nR, trial, npos, nzero
    real colvector  xx, /*ww,*/ l, r, L, R, xpos, xneg

    if (rows(x)<=1) return(0)      // returning zero if n=1
    xpos    = select(x, x:>0)      // obervations > median
    npos    = rows(xpos)           // number of obs > median
    xneg    = select(x, x:<0)      // observations < median
    nzero   = sum(x:==0)           // number of obs = median
    n       = npos + nzero         // number of rows in search matrix
    q       = rows(xneg) + nzero   // number of columns in search matrix
    xx      = /*ww =*/ J(n, 1, .)  // temp vector for matrix elements
    l = L   = J(n, 1, 1)           // indices of left boundary (old and new)
    r = R   = J(n, 1, q)           // indices of right boundary (old and new)
    nl = nL = 0                    // number of cells below left boundary
    nr = nR = n * q                // number of cells within right boundary
    k       = n*q/2                // target quantile
    while ((nr-nl)>n) {
        // get trial value
        m = 0
        for (i=1; i<=n; i++) {
            if (l[i]<=r[i]) {
                // high median within row
                xx[++m] = -__mm_mc_el(xpos, xneg, npos, nzero, i,
                           l[i]+trunc((r[i]-l[i]+1)/2))
                /*m++
                ww[m] = r[i] - l[i] + 1
                xx[m] = -__mm_mc_el(xpos, xneg, npos, nzero, i,
                          l[i]+trunc((ww[m])/2))*/
            }
        }
        trial = _mm_hl_qhi(xx[|1 \ m|], .5)
        /*the unweighted quantile is faster; results are the same*/
        /*trial = _mm_hl_qhi_w(xx[|1 \ m|], ww[|1 \ m|], .5)*/
        // move right border
        j = 0
        for (i=n; i>=1; i--) {
            if (j<q) {
                while (-__mm_mc_el(xpos, xneg, npos, nzero, i, j+1)<trial) {
                    j++
                    if (j==q) break
                }
            }
            R[i] = j
        }
        nR = sum(R)
        if (nR>k) {
            swap(r, R)
            nr = nR
            continue
        }
        // move left border
        j = q + 1
        for (i=1; i<=n; i++) {
            while (-__mm_mc_el(xpos, xneg, npos, nzero, i, j-1)>trial) {
                j--
            }
            L[i] = j
        }
        nL = sum(L) - n
        if (nL<k) {
            swap(l, L)
            nl = nL
            continue
        }
        // trial = low quantile = high quantile
        if (ceil(k)!=k | (nR<k & nL>k)) return(-trial)
        // trial = low quantile
        if (nL==k) {
            m = 0
            for (i=1; i<=n; i++) {
                if (L[i]>q) continue
                xx[++m] = -__mm_mc_el(xpos, xneg, npos, nzero, i, L[i])
            }
            return(-(trial+min(xx[|1 \ m|]))/2)
        }
        // trial = high quantile
        for (i=1; i<=n; i++) {
            xx[i] = -__mm_mc_el(xpos, xneg, npos, nzero, i, R[i])
        }
        return(-(trial+max(xx))/2)
    }
    // get target value from remaining candidates
    m = 0
    for (i=1; i<=n; i++) {
        if (l[i]<=r[i]) {
            for (j=l[i]; j<=r[i]; j++) {
                m++
                xx[m] = -__mm_mc_el(xpos, xneg, npos, nzero, i, j)
            }
        }
    }
    return(-_mm_hl_q(xx[|1 \ m|], k, nl))
}

real scalar __mm_mc_el(real colvector xpos, real colvector xneg, 
    real scalar npos, real scalar nzero, real scalar i, real scalar j)
{
    if (i<=npos) {
        if (j<=nzero) return(1)
        // => (j>nzero)
        return((xpos[i] + xneg[j-nzero])/(xpos[i] - xneg[j-nzero]))
    }
    // => (i>npos)
    if (j>nzero)  return(-1)
    // => (j<=nzero)
    return(sign((npos+nzero-i+1)-j))
}

real scalar _mm_mc_w(real colvector x, real colvector w) // assumes med(x)=0
{
    real scalar     i, j, m, k, n, q, nl, nr, trial, npos, nzero, Wl, WR, WL, W
    real colvector  xx, ww, l, r, L, R, p, ccw, xpos, xneg, wpos, wneg, wzero

    if (rows(x)<=1) return(0)      // returning zero if n=1
    p = (x:>0)
    xpos    = select(x, p)         // obervations > median
    wpos    = select(w, p)         // weights of obervations >= median
    npos    = rows(xpos)           // number of obs > median
    p = (x:<0)
    xneg    = select(x, p)         // observations < median
    wneg    = select(w, p)         // weights of obervations <= median
    wzero   = select(w, x:==0)     // weights of obervations = median
    nzero   = rows(wzero)          // number of obs = median
    if (nzero>0) {
        wpos = wpos \ wzero
        wneg = wzero[nzero::1] \ wneg // need to use reverse ordered wzero
    }
    n       = npos + nzero         // number of rows in search matrix
    q       = rows(xneg) + nzero   // number of columns in search matrix
    xx = ww = J(n, 1, .)           // temp vector for matrix elements
    l = L   = J(n, 1, 1)           // indices of left boundary (old and new)
    r = R   = J(n, 1, q)           // indices of right boundary (old and new)
    nl      = 0                    // number of cells below left boundary
    nr      = n * q                // number of cells within right boundary
    ccw     = quadrunningsum(wneg) // cumulative column weights
    W       = quadsum(wpos:*ccw[rows(ccw)]) // sum of weights in search matrix
    Wl = WL = 0                    // sum of weights below left boundary
    WR      = W                    // sum of weights within right boundary
    k       = W/2                  // target quantile
    while ((nr-nl)>n) {
        // get trial value
        m = 0
        for (i=1; i<=n; i++) {
            if (l[i]<=r[i]) {
                // high median within row
                xx[++m] = -__mm_mc_el(xpos, xneg, npos, nzero, i,
                           l[i]+trunc((r[i]-l[i]+1)/2))
            }
        }
        trial = _mm_hl_qhi(xx[|1 \ m|], .5)
        // move right border
        j = 0
        for (i=n; i>=1; i--) {
            if (j<q) {
                while (-__mm_mc_el(xpos, xneg, npos, nzero, i, j+1)<trial) {
                    j++
                    if (j==q) break
                }
            }
            R[i] = j
        }
        p = (R:>0)
        if (any(p)) WR = quadsum(select(wpos, p) :* ccw[select(R, p)])
        else        WR = 0
        if (WR>k) {
            swap(r, R)
            nr = sum(R)
            continue
        }
        // move left border
        j = q + 1
        for (i=1; i<=n; i++) {
            while (-__mm_mc_el(xpos, xneg, npos, nzero, i, j-1)>trial) {
                j--
            }
            L[i] = j
        }
        p = (L:>1)
        WL = quadsum(select(wpos, p) :* ccw[select(L, p):-1])
        if (WL<k) {
            swap(l, L)
            Wl = WL
            nl = sum(L) - n
            continue
        }
        // trial = low quantile = high quantile
        if (WR==WL | (WR<k & WL>k)) return(-trial)
        // trial = low quantile
        if (WL==k) {
            m = 0
            for (i=1; i<=n; i++) {
                if (L[i]>q) continue
                xx[++m] = -__mm_mc_el(xpos, xneg, npos, nzero, i, L[i])
            }
            return(-(trial+min(xx[|1 \ m|]))/2)
        }
        // trial = high quantile
        for (i=1; i<=n; i++) {
            xx[i] = -__mm_mc_el(xpos, xneg, npos, nzero, i, R[i])
        }
        return(-(trial+max(xx))/2)
    }
    // get target value from remaining candidates
    m = 0
    for (i=1; i<=n; i++) {
        if (l[i]<=r[i]) {
            for (j=l[i]; j<=r[i]; j++) {
                m++
                xx[m] = -__mm_mc_el(xpos, xneg, npos, nzero, i, j)
                ww[m] = wpos[i] * wneg[j]
            }
        }
    }
    return(-_mm_hl_q_w(xx[|1 \ m|], ww[|1 \ m|], k, Wl))
}

real scalar _mm_mc_fw(real colvector x, real colvector w) // assumes med(x)=0
{
    real scalar     i, j, m, k, n, q, nl, nr, trial, npos, nzero, Wl, WR, WL, W
    real colvector  xx, ww, l, r, L, R, p, ccw, ccwz, xpos, xneg, wpos, wneg, wnegz, wzero

    if (any(trunc(w):!=w)) _error(3498, "non-integer frequency not allowed")
    if (rows(x)<=1) return(0)       // returning zero if n=1
    p = (x:>0)
    xpos    = select(x, p)          // obervations > median
    wpos    = select(w, p)          // weights of obervations >= median
    npos    = rows(xpos)            // number of obs > median
    p = (x:<0)
    xneg    = select(x, p)          // observations < median
    wneg    = select(w, p)          // weights of obervations <= median
    wzero   = sum(select(w, x:==0)) // aggregate weights for x==median
    nzero   = (wzero>0)             // has x==median
    if (nzero>0) {
        wpos  = wpos \ 1
        wnegz = (wzero>1 ? (comb(wzero, 2), wzero, comb(wzero, 2))' : wzero) \ wneg*wzero
        wneg  = wzero \ wneg
    }
    n       = npos + nzero          // number of rows in search matrix
    q       = nzero + rows(xneg)    // number of columns in search matrix
    xx = ww = J(n, 1, .)            // temp vector for matrix elements
    l = L   = J(n, 1, 1)            // indices of left boundary (old and new)
    r = R   = J(n, 1, q)            // indices of right boundary (old and new)
    if (wzero>1) {
        r[n] = q+2
        R[n] = q+2
    }
    nl      = 0                     // number of cells below left boundary
    nr      = sum(R)                // number of cells within right boundary
    ccw     = runningsum(wneg)      // cumulative column weights
    if (nzero==0) {
        W   = quadsum(wpos:*ccw[rows(ccw)]) // sum of weights in search matrix
    }
    else {
        ccwz = runningsum(wnegz)    // cumulative column weights in last row
        if (npos>=1) W = quadsum(wpos[|1 \ n-1|]:*ccw[rows(ccw)]) + ccwz[rows(ccwz)]
        else         W = ccwz[rows(ccwz)]
    }
    Wl = WL = 0                       // sum of weights below left boundary
    WR      = W                       // sum of weights within right boundary
    k       = W/2                     // target quantile
    while ((nr-nl)>n) {
        // get trial value
        m = 0
        for (i=1; i<=(n-(wzero>1)); i++) {
            if (l[i]<=r[i]) {
                // high median within row
                xx[++m] = -__mm_mc_el(xpos, xneg, npos, nzero, i,
                           l[i]+trunc((r[i]-l[i]+1)/2))
            }
        }
        if (wzero>1) { // handle last row
            if (l[n]<=r[n]) {
                m++
                xx[m] = -__mm_mc_el(xpos, xneg, npos-1, 3, n,
                           l[n]+trunc((r[n]-l[n]+1)/2))
            }
        }
        trial = _mm_hl_qhi(xx[|1 \ m|], .5)
        // move right border
        if (wzero>1) { // handle last row
            j = 0
            while (-__mm_mc_el(xpos, xneg, npos-1, 3, n, j+1)<trial) {
                j++
                if (j==(q+2)) break
            }
            R[n] = j
        }
        j = 0
        for (i=(n-(wzero>1)); i>=1; i--) {
            if (j<q) {
                while (-__mm_mc_el(xpos, xneg, npos, nzero, i, j+1)<trial) {
                    j++
                    if (j==q) break
                }
            }
            R[i] = j
        }
        p = (R:>0)
        if (nzero>0) p[n] = 0
        if (any(p)) WR = quadsum(select(wpos, p) :* ccw[select(R, p)])
        else        WR = 0
        if (nzero>0) {
            if (R[n]>0) WR = WR + ccwz[R[n]]
        }
        if (WR>k) {
            swap(r, R)
            nr = sum(R)
            continue
        }
        // move left border
        j = q + 1
        for (i=1; i<=(n-(wzero>1)); i++) {
            while (-__mm_mc_el(xpos, xneg, npos, nzero, i, j-1)>trial) {
                j--
            }
            L[i] = j
        }
        if (wzero>1) { // handle last row
            j = q + 3
            while (-__mm_mc_el(xpos, xneg, npos-1, 3, n, j-1)>trial) {
                j--
            }
            L[n] = j
        }
        p = (L:>1)
        if (nzero>0) p[n] = 0
        if (any(p)) WL = quadsum(select(wpos, p) :* 
            (rows(ccw)==1 ? ccw[select(L, p):-1]' : ccw[select(L, p):-1]))
        else        WL = 0
        if (nzero>0) {
            if (L[n]>1) WL = WL + ccwz[L[n]-1]
        }
        if (WL<k) {
            swap(l, L)
            Wl = WL
            nl = sum(L) - n
            continue
        }
        // trial = low quantile = high quantile
        if (ceil(k)!=k | (WR<k & WL>k)) return(-trial)
        // trial = low quantile
        if (WL==k) {
            m = 0
            for (i=1; i<=(n-(wzero>1)); i++) {
                if (L[i]>q) continue
                xx[++m] = -__mm_mc_el(xpos, xneg, npos, nzero, i, L[i])
            }
            if (wzero>1) { // handle last row
                if (L[n]<=(q+2)) {
                    xx[++m] = -__mm_mc_el(xpos, xneg, npos-1, 3, n, L[n])
                }
            }
            return(-(trial+min(xx[|1 \ m|]))/2)
        }
        // trial = high quantile
        for (i=1; i<=(n-(wzero>1)); i++) {
            xx[i] = -__mm_mc_el(xpos, xneg, npos, nzero, i, R[i])
        }
        if (wzero>1) { // handle last row
            xx[n] = -__mm_mc_el(xpos, xneg, npos-1, 3, n, R[n])
        }
        return(-(trial+max(xx))/2)
    }
    // get target value from remaining candidates
    m = 0
    for (i=1; i<=(n-nzero); i++) {
        if (l[i]<=r[i]) {
            for (j=l[i]; j<=r[i]; j++) {
                m++
                xx[m] = -__mm_mc_el(xpos, xneg, npos, nzero, i, j)
                ww[m] = wpos[i] * wneg[j]
            }
        }
    }
    if (nzero>0) { // handle last row
        if (l[n]<=r[n]) {
            for (j=l[n]; j<=r[n]; j++) {
                m++
                xx[m] = -__mm_mc_el(xpos, xneg, npos-(wzero>1),
                    (wzero>1 ? 3 : 1), n, j)
                ww[m] = wnegz[j]
            }
        }
    }
    return(-_mm_hl_q_w(xx[|1 \ m|], ww[|1 \ m|], k, Wl))
}

real scalar _mm_mc_naive(real colvector x) // noweights; assumes med(x)=0
{
    real scalar    i, j, m, n, q, npos, nzero
    real colvector xx, xpos, xneg 

    if (rows(x)<=1) return(0) // returning zero if n=1
    xpos   = select(x, x:>0)
    npos   = rows(xpos)
    xneg   = select(x, x:<0)
    nzero  = sum(x:==0)
    n      = npos + nzero
    q      = rows(xneg) + nzero
    m = 0
    xx = J(n*q, 1, .)
    for (i=1; i<=n; i++) {
        for (j=1; j<=q; j++) {
            xx[++m] = __mm_mc_el(xpos, xneg, npos, nzero, i, j)
        }
    }
    return(_mm_hl_q(xx, .5))
}

real scalar _mm_mc_naive_w(real colvector x, real colvector w) // assumes med(x)=0
{
    real scalar    i, j, m, n, q, npos, nzero
    real colvector xx, ww, xpos, xneg, wpos, wneg, wzero

    if (rows(x)<=1) return(0) // returning zero if n=1
    xpos   = select(x, x:>0)
    wpos   = select(w, x:>0)
    npos   = rows(xpos)
    xneg   = select(x, x:<0)
    wneg   = select(w, x:<0)
    wzero  = select(w, x:==0)
    nzero  = rows(wzero)
    if (nzero>0) {
        wpos = wpos \ wzero
        wneg = wzero[nzero::1] \ wneg // need to use reverse ordered wzero
    }
    n      = npos + nzero
    q      = rows(xneg) + nzero
    m = 0
    xx = ww = J(n*q, 1, .)
    for (i=1; i<=n; i++) {
        for (j=1; j<=q; j++) {
            m++
            xx[m] = __mm_mc_el(xpos, xneg, npos, nzero, i, j)
            ww[m] = wpos[i] * wneg[j]
        }
    }
    return(_mm_hl_q_w(xx, ww, .5))
}

real scalar _mm_mc_naive_fw(real colvector x, real colvector w) // assumes med(x)=0
{
    real scalar    i, j, m, n, q, npos, nzero, wzero
    real colvector xx, ww, xpos, xneg, wpos, wneg

    if (any(trunc(w):!=w)) _error(3498, "non-integer frequency not allowed")
    if (rows(x)<=1) return(0) // returning zero if n=1
    xpos   = select(x, x:>0)
    wpos   = select(w, x:>0)
    npos   = rows(xpos)
    xneg   = select(x, x:<0)
    wneg   = select(w, x:<0)
    wzero  = sum(select(w, x:==0)) // aggregate weights for x==median
    nzero  = (wzero>0)
    if (nzero>0) {
        wpos = wpos \ wzero
        wneg = wzero \ wneg
    }
    n      = npos + nzero
    q      = rows(xneg) + nzero
    m = 0
    xx = ww = J(n*q + 2*(wzero>1), 1, .)
    for (i=1; i<=n; i++) {
        for (j=1; j<=q; j++) {
            m++
            if (i>npos & j==nzero) { // x==median
                xx[m] = 0
                ww[m] = wzero
                if (wzero>1) {
                    m++
                    xx[m] = 1
                    ww[m] = comb(wzero, 2)
                    m++
                    xx[m] = -1
                    ww[m] = comb(wzero, 2)
                }
            }
            else {
                xx[m] = __mm_mc_el(xpos, xneg, npos, nzero, i, j)
                ww[m] = wpos[i] * wneg[j]
            }
        }
    }
    return(_mm_hl_q_w(xx, ww, .5))
}

// helper functions for quantiles

real scalar _mm_hl_q(real colvector x, real scalar P, 
    | real scalar offset) // must be integer; changes meaning of P if specified
{   // quantile (definition 2)
    real scalar    j0, j1, n, k
    real colvector p

    n = rows(x)
    if (n<1) return(.)
    if (n==1) return(x)
    if (args()==3) k = P     // P is a count (possibly noninteger)
    else           k = P * n // P is a proportion
    j0 = ceil(k)      - (args()==3 ? offset : 0) // index of low quantile
    j1 = floor(k) + 1 - (args()==3 ? offset : 0) // index of high quantile
    if (j0<1)      j0 = 1
    else if (j0>n) j0 = n
    if (j1<1)      j1 = 1
    else if (j1>n) j1 = n
    p = order(x, 1)
    if (j0==j1) return(x[p[j1]])
    return((x[p[j0]] + x[p[j1]])/2)
}

real scalar _mm_hl_qhi(real colvector x, real scalar P)
{   // high quantile
    real scalar    j, n
    real colvector p

    n = rows(x)
    if (n<1) return(.)
    p = order(x, 1)
    j = floor(P * n) + 1
    if (j<1)      j = 1
    else if (j>n) j = n
    return(x[p[j]])
}

real scalar _mm_hl_q_w(real colvector x, real colvector w, real scalar P, 
    | real scalar offset) // changes meaning of P if specified
{
    real scalar    n, i, k
    real colvector p, cw

    n = rows(x)
    if (n<1) return(.)
    p = order(x, 1)
    if (anyof(w, 0)) {
         p = select(p, w[p]:!=0)
         n = rows(p)
    }
    if (n<1) return(.)
    if (n==1) return(x[p])
    if (args()==4) {
        cw = quadrunningsum(offset \ w[p])[|2 \ n+1|]
        k = P // P is a count
    }
    else {
        cw = quadrunningsum(w[p])
        k = P * cw[n] // P is a proportion
    }
    if (k>=cw[n]) return(x[p[n]])
    for (i=1; i<=n; i++) {
        if (k>cw[i]) continue
        if (k==cw[i]) return((x[p[i]]+x[p[i+1]])/2)
        return(x[p[i]])
    }
    // cannot be reached
}

/*
real scalar _mm_hl_qhi_w(real colvector x, real colvector w, real scalar P)
{   // high quantile (weighted)
    real scalar    i, n, k
    real colvector p, cw

    n = rows(x)
    if (n<1) return(.)
    p = order(x, 1)
    cw = quadrunningsum(w[p])
    k  = cw[n] * P
    if (k>=cw[n]) return(x[p[n]])
    for (i=1; i<=n; i++) {
        if (k>=cw[i]) continue
        return(x[p[i]])
    }
    // cannot be reached
}
*/

end
