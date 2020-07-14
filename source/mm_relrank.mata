*! version 1.1.0  13jul2020  Ben Jann
version 9.2
mata:

real matrix mm_relrank(
    real matrix X,
    real colvector w,
    real matrix Y,
    | real scalar mid,
      real scalar nonorm, 
      real scalar brk,
      real colvector w2)
{
    if (args()<4) mid = 0
    if (args()<5) nonorm = 0
    if (args()<6) brk = 0
    if (args()<7) w2 = 1
    if (rows(w)!=1 & rows(w)!=rows(X)) _error(3200)
    if (rows(w2)!=1 & rows(w2)!=rows(Y)) _error(3200)
    if (cols(X)==1 & cols(Y)!=1 & rows(Y)==1 & rows(w2)==1)
        return(_mm_relrank_sort(X, w, Y', mid, nonorm, brk, w2)')
    return(_mm_relrank_sort(X, w, Y, mid, nonorm, brk, w2))
}

real matrix _mm_relrank_sort(
    real matrix X,
    real colvector w,
    real matrix Y,
    real scalar mid,
    real scalar nonorm, 
    real scalar brk,
    real colvector w2)
{
    real scalar    i, c, c1, c2
    real colvector p, sX, sw, p2, sY, sw2
    real matrix    R

    c1 = cols(X); c2 = cols(Y)
    c = max((c1,c2))
    R = J(rows(Y), c, .)
    if (c1==c2) {
        for (i=c; i; i--) {
            if (rows(w)==1) {; p = order(X[,i],1); sX = X[p,i]; sw = w; }
            else {; p = order((X[,i],w),(1,2)); sX = X[p,i]; sw = w[p]; }
            if (rows(w2)==1) {; p2 = order(Y[,i],1); sY = Y[p2,i]; sw2 = w2; }
            else {; p2 = order((Y[,i],w2),(1,2)); sY = Y[p2,i]; sw2 = w2[p2]; }
            R[p2,i] = _mm_relrank(sX, sw, sY, mid, nonorm, brk, sw2)
        }
        return(R)
    }
    if (c1==1) {
        if (rows(w)==1) {; p = order(X,1); sX = X[p]; sw = w; }
        else {; p = order((X,w),(1,2)); sX = X[p]; sw = w[p]; }
        for (i=c; i; i--) {
            if (rows(w2)==1) {; p2 = order(Y[,i],1); sY = Y[p2,i]; sw2 = w2; }
            else {; p2 = order((Y[,i],w2),(1,2)); sY = Y[p2,i]; sw2 = w2[p2]; }
            R[p2,i] = _mm_relrank(sX, sw, sY, mid, nonorm, brk, sw2)
        }
        return(R)
    }
    if (c2==1) {
        if (rows(w2)==1) {; p2 = order(Y,1); sY = Y[p2]; sw2 = w2; }
        else {; p2 = order((Y,w2),(1,2)); sY = Y[p2]; sw2 = w2[p2]; }
        for (i=c; i; i--) {
            if (rows(w)==1) {; p = order(X[,i],1); sX = X[p,i]; sw = w; }
            else {; p = order((X[,i],w),(1,2)); sX = X[p,i]; sw = w[p]; }
            R[p2,i] = _mm_relrank(sX, sw, sY, mid, nonorm, brk, sw2)
        }
        return(R)
    }
    _error(3200)
}

real colvector _mm_relrank(
    real colvector X,
    real colvector w,
    real colvector Y,
    | real scalar mid,
      real scalar nonorm, 
      real scalar brk,
      real colvector w2)
{
    real matrix cdf
    
    if (args()<4) mid = 0
    if (args()<5) nonorm = 0
    if (args()<6) brk = 0
    if (args()<7) w2 = 1
    if (rows(X)==0 | rows(Y)==0) return(J(rows(Y),1,.))
    
    // obtain CDF at unique values of X
    cdf = _mm_ecdf2(X, w, 0, nonorm)
    
    // compute relative ranks
    if (brk==0) {
        if (mid==0) return(_mm_relrank_1(Y, cdf[,1], cdf[,2]))
                    return(_mm_relrank_2(Y, cdf[,1], cdf[,2]))
    }
    if (rows(w2)==1) {
        if (mid==0) return(_mm_relrank_3(Y, cdf[,1], cdf[,2]))
                    return(_mm_relrank_4(Y, cdf[,1], cdf[,2]))
    }
    if (mid==0)     return(_mm_relrank_5(Y, w2, cdf[,1], cdf[,2]))
                    return(_mm_relrank_6(Y, w2, cdf[,1], cdf[,2]))
}

real colvector _mm_relrank_1(real colvector x, real colvector y, real colvector cdf)
{    // case 1: brk = 0, mid = 0
    real scalar    i, j, xi
    real colvector r
    
    i = rows(x)
    r = J(i, 1, 0)
    j = rows(y)
    for (; i; i--) {
        xi = x[i]
        for (; j; j--) {
            if (y[j]<=xi) break
        }
        if (j) r[i] = cdf[j]
        else break // x[i] is smaller than min(y)
    }
    return(r)
}
real colvector _mm_relrank_2(real colvector x, real colvector y, real colvector cdf)
{   // case 2: brk = 0, mid = 1
    real scalar    i, j, xi
    real colvector r, step
    
    i = rows(x)
    r = J(i, 1, 0)
    j = rows(y)
    step = mm_diff(0\cdf)
    for (; i; i--) {
        xi = x[i]
        for (; j; j--) {
            if (y[j]<=xi) break
        }
        if (j) {
            if (y[j]==xi) r[i] = cdf[j] - step[j]/2
            else          r[i] = cdf[j]
        }
        else break // x[i] is smaller than min(y)
    }
    return(r)
}
real colvector _mm_relrank_3(real colvector x, real colvector y, real colvector cdf)
{
    // case 3: brk = 1, mid = 0, no weights
    real scalar    i, j, k, xi
    real colvector r, step

    i = rows(x)
    r = J(i, 1, 0)
    j = rows(y)
    step = mm_diff(0\cdf)
    for (; i; i--) {
        xi = x[i]
        for (; j; j--) {
            if (y[j]<=xi) break
        }
        if (j) {
            r[i] = cdf[j]
            if (y[j]==xi) {
                for (k=i-1; k; k--) { // find ties in x
                    if (x[k]<xi) break
                }
                if ((++k)==i) continue // no ties
                r[|k\i-1|] = cdf[j] :- step[j] :* (i:-(k::i-1)) / (i-k+1)
                i = k
            }
        }
        else break // x[i] is smaller than min(y)
    }
    return(r)
}
real colvector _mm_relrank_4(real colvector x, real colvector y, real colvector cdf)
{
    // case 4: brk = 1, mid = 1, no weights
    real scalar    i, j, k, xi
    real colvector r, step
    
    i = rows(x)
    r = J(i, 1, 0)
    j = rows(y)
    step = mm_diff(0\cdf)
    for (; i; i--) {
        xi = x[i]
        for (; j; j--) {
            if (y[j]<=xi) break
        }
        if (j) {
            if (y[j]==xi) {
                for (k=i-1; k; k--) { // find ties in x
                    if (x[k]<xi) break
                }
                if ((++k)==i) {
                    r[i] = cdf[j] - step[j] * 0.5
                    continue
                }
                r[|k\i|] = cdf[j] :- step[j] :* ((i+.5):-(k::i)) / (i-k+1)
                i = k
            }
            else r[i] = cdf[j]
        }
        else break // x[i] is smaller than min(y)
    }
    return(r)
}
real colvector _mm_relrank_5(real colvector x, real colvector w, 
    real colvector y, real colvector cdf)
{
    // case 5: brk = 1, mid = 0, weighted
    real scalar    i, j, k, xi, W
    real colvector r, step, ww
    
    i = rows(x)
    r = J(i, 1, 0)
    j = rows(y)
    step = mm_diff(0\cdf)
    for (; i; i--) {
        xi = x[i]
        for (; j; j--) {
            if (y[j]<=xi) break
        }
        if (j) {
            r[i] = cdf[j]
            if (y[j]==xi) {
                for (k=i-1; k; k--) { // find ties in x
                    if (x[k]<xi) break
                }
                if ((++k)==i) continue // no ties
                ww = mm_colrunsum(w[|k\i|])
                W  = ww[rows(ww)]
                if (W==0) {
                    r[|k\i-1|] = cdf[j] :- step[j] :* (i:-(k::i-1)) / (i-k+1)
                }
                else {
                    ww = ww[|1\rows(ww)-1|]
                    r[|k\i-1|] = cdf[j] :- step[j] :* (W:-ww) / W
                }
                i = k
            }
        }
        else break // x[i] is smaller than min(y)
    }
    return(r)
}
real colvector _mm_relrank_6(real colvector x, real colvector w, 
    real colvector y, real colvector cdf)
{
    // case 6: brk = 1, mid = 1, weighted
    real scalar    i, j, k, xi, W
    real colvector r, step, ww
    
    i = rows(x)
    r = J(i, 1, 0)
    j = rows(y)
    step = mm_diff(0\cdf)
    for (; i; i--) {
        xi = x[i]
        for (; j; j--) {
            if (y[j]<=xi) break
        }
        if (j) {
            if (y[j]==xi) {
                for (k=i-1; k; k--) { // find ties in x
                    if (x[k]<xi) break
                }
                if ((++k)==i) {
                    r[i] = cdf[j] - step[j] * 0.5
                    continue // no ties
                }
                ww = mm_colrunsum(w[|k\i|])
                W  = ww[rows(ww)]
                if (W==0) {
                    // if all observations have zero weight
                    r[|k\i|] = cdf[j] :- step[j] :* ((i+.5):-(k::i)) / (i-k+1)
                }
                else {
                    ww = ww - 0.5 * w[|k\i|]
                    r[|k\i|] = cdf[j] :- step[j] :* (W:-ww) / W
                }
                i = k
            }
            else r[i] = cdf[j]
        }
        else break // x[i] is smaller than min(y)
    }
    return(r)
}

end
