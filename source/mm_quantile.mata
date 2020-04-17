*! version 1.0.8  20dec2007  Ben Jann
version 9.2
mata:

real matrix mm_quantile(real matrix X, | real colvector w,
 real matrix P, real scalar altdef)
{
    real rowvector result
    real scalar c, cX, cP, r, i

    if (args()<2) w = 1
    if (args()<3) P = (0, .25, .50, .75, 1)'
    if (args()<4) altdef = 0
    if (cols(X)==1 & cols(P)!=1 & rows(P)==1)
     return(mm_quantile(X, w, P', altdef)')
    if (missing(P) | missing(X) | missing(w)) _error(3351)
    if (rows(w)!=1 & rows(w)!=rows(X)) _error(3200)
    r = rows(P)
    c = max(((cX=cols(X)), (cP=cols(P))))
    if (cX!=1 & cX<c) _error(3200)
    if (cP!=1 & cP<c) _error(3200)
    if (rows(X)==0 | r==0 | c==0) return(J(r,c,.))
    if (c==1) return(_mm_quantile(X, w, P, altdef))
    result = J(r, c, .)
    if (cP==1) for (i=1; i<=c; i++)
     result[,i] = _mm_quantile(X[,i], w, P, altdef)
    else if (cX==1) for (i=1; i<=c; i++)
     result[,i] = _mm_quantile(X, w, P[,i], altdef)
    else for (i=1; i<=c; i++)
     result[,i] = _mm_quantile(X[,i], w, P[,i], altdef)
    return(result)
}

real colvector _mm_quantile(
 real colvector X,
 real colvector w,
 real colvector P,
 real scalar altdef)
{
    real colvector g, j, j1, p
    real scalar N

    if (w!=1) return(_mm_quantilew(X, w, P, altdef))
    N = rows(X)
    p = order(X,1)
    if (altdef) g = P*N + P
    else g = P*N
    j = floor(g)
    if (altdef) g = g - j
    else g = 0.5 :+ 0.5*((g - j):>0)
    j1 = j:+1
    j = j :* (j:>=1)
    _editvalue(j, 0, 1)
    j = j :* (j:<=N)
    _editvalue(j, 0, N)
    j1 = j1 :* (j1:>=1)
    _editvalue(j1, 0, 1)
    j1 = j1 :* (j1:<=N)
    _editvalue(j1, 0, N)
    return((1:-g):*X[p[j]] + g:*X[p[j1]])
}

real colvector _mm_quantilew(
 real colvector X,
 real colvector w,
 real colvector P,
 real scalar altdef)
{
    real colvector Q, pi, pj
    real scalar i, I, j, jj, J, rsum, W
    pointer scalar ww

    I  = rows(X)
    ww = (rows(w)==1 ? &J(I,1,w) : &w)
    if (altdef) return(_mm_quantilewalt(X, *ww, P))
    W  = quadsum(*ww)
    pi = order(X, 1)
    if (anyof(*ww, 0)) {
        pi = select(pi,(*ww)[pi]:!=0)
        I = rows(pi)
    }
    pj = order(P, 1)
    J  = rows(P)
    Q  = J(J, 1, .)
    j  = 1
    jj = pj[1]
    rsum = 0
    for (i=1; i<=I; i++) {
        rsum = rsum + (*ww)[pi[i]]
        if (i<I) {
            if (rsum<P[jj]*W) continue
            if (X[pi[i]]==X[pi[i+1]]) continue
        }
        while (1) {
            if (rsum>P[jj]*W | i==I) Q[jj] = X[pi[i]]
            else Q[jj] = (X[pi[i]] + X[pi[i+1]])/2
            j++
            if (j>J) break
            jj = pj[j]
            if (i<I & rsum<P[jj]*W) break
        }
        if (j>J) break
    }
    return(Q)
}

real colvector _mm_quantilewalt(
 real colvector X,
 real colvector w,
 real colvector P)
{
    real colvector Q, pi, pj
    real scalar i, I, j, jj, J, rsum, rsum0, W, ub, g

    W  = quadsum(w) + 1
    pi = order(X, 1)
    if (anyof(w, 0)) pi = select(pi, w[pi]:!=0)
    I  = rows(pi)
    pj = order(P, 1)
    J  = rows(P)
    Q  = J(J, 1, .)
    rsum = w[pi[1]]
    for (j=1; j<=J; j++) {
        jj = pj[j]
        if (P[jj]*W <= rsum) Q[jj] = X[pi[1]]
        else break
    }
    for (i=2; i<=I; i++) {
        rsum0 = rsum
        rsum = rsum + w[pi[i]]
        if (i<I & rsum < P[jj]*W) continue
        while (1) {
            ub = rsum0+1
            if (P[jj]*W>=ub | X[pi[i]]==X[pi[i-1]]) Q[jj] = X[pi[i]]
            else {
                g = (ub - P[jj]*W) / (ub - rsum0)
                Q[jj] = X[pi[i-1]]*g + X[pi[i]]*(1-g)
            }
            j++
            if (j>J) break
            jj = pj[j]
            if (i<I & rsum < P[jj]*W) break
        }
        if (j>J) break
    }
    return(Q)
}
end
