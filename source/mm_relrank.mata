*! version 1.0.7  03aug2007  Ben Jann
version 9.2
mata:

real matrix mm_relrank(
    real matrix X,
    real colvector w,
    real matrix Q,
    | real scalar mid)
{
    real rowvector result
    real scalar c, cX, cQ, r, i

    if (args()<4) mid = 0

    if (cols(X)==1 & cols(Q)!=1 & rows(Q)==1) return(mm_relrank(X, w, Q', mid)')
    if (rows(w)!=1 & rows(w)!=rows(X)) _error(3200)
    r = rows(Q)
    c = max(((cX=cols(X)), (cQ=cols(Q))))
    if (cX!=1 & cX<c) _error(3200)
    if (cQ!=1 & cQ<c) _error(3200)
    if (rows(X)==0 | r==0 | c==0) return(J(r,c,.))
    if (c==1) return(_mm_relrank(X, w, Q, mid))
    result = J(r, c, .)
    if (cQ==1)
     for (i=1; i<=c; i++) result[,i] = _mm_relrank(X[,i], w, Q, mid)
    else if (cX==1)
     for (i=1; i<=c; i++) result[,i] = _mm_relrank(X, w, Q[,i], mid)
    else
     for (i=1; i<=c; i++) result[,i] = _mm_relrank(X[,i], w, Q[,i], mid)
    return(result)
}

real colvector _mm_relrank(
 real colvector X,
 real colvector w,
 real colvector Q,
 real scalar mid)
{
    real colvector P, cdf, pi, pj, sortx
    real scalar i, I, j, J

    I = rows(X)
    pi = order((X,(1::I)),(1,2)) //stable sort order
    sortx = X[pi]
    cdf = __mm_ranks(sortx, (rows(w)!=1 ? w[pi] : w), 3, 0, 1)

    J = rows(Q)
    pj = order(Q, 1)
    P = J(J, 1, 0)
    i = I
    for(j=J; j>=1; j--) {
        for (; i>=1; i--) {
            if (sortx[i]<=Q[pj[j]]) break
        }
        if (i<1) break
        P[pj[j]] = cdf[i]
    }
    if (mid==0) return(P)

    cdf = _mm_ranks_mid(sortx, cdf, mid)
    i = 1
    for (j=1; j<=J; j++) {
        for (; i<I; i++) {
            if (sortx[i]>=Q[pj[j]]) break
        }
        if (sortx[i]==Q[pj[j]]) P[pj[j]] = cdf[i]
    }
    return(P)
}
end
