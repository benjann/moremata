*! version 1.1.0  19oct2020  Ben Jann
version 9.2
mata:

real matrix mm_ranks(real matrix X, | real colvector w,
 real scalar method, real scalar mid, real scalar norm)
{
    real scalar c, r, i
    real matrix result

    if (args()<2) w = 1
    if (args()<3) method = 0
    if (args()<4) mid = 0
    if (args()<5) norm = 0
    c = cols(X)
    r = rows(X)
    if (rows(w)!=1 & rows(w)!=r) _error(3200)
    if (c<1) return(J(r,0,.))
    if (c==1) return(_mm_ranks_sort(X, w, method, mid, norm))
    result = J(r, c, .)
    for (i=1; i<=c; i++) {
        result[,i] = _mm_ranks_sort(X[,i], w, method, mid, norm)
    }
    return(result)
}

real colvector _mm_ranks_sort(real colvector x, real colvector w,
 real scalar method, real scalar mid, real scalar norm)
{
    real scalar    I
    real colvector p, r

    I = rows(x)
    if (I==0) return(J(0,1,.))
    if      (method==4 & rows(w)!=1) p = order((x,w),(1,2))
    else if (method==0 | method==4)  p = order(x,1)
    else if (rows(w)!=1) p = order((x,w),(1,2)) // include w for stable results
    else p = order(x,1)
    r = __mm_ranks(x[p], (rows(w)!=1 ? w[p] : w), method, mid, norm)
    r[p] = r
    return(r)
}

real colvector _mm_ranks(real colvector x, | real colvector w,
 real scalar method, real scalar mid, real scalar norm) // sorted input assumed
{
    if (args()<2) w = 1
    if (args()<3) method = 0
    if (args()<4) mid = 0
    if (args()<5) norm = 0
    if (rows(w)!=1 & rows(w)!=rows(x)) _error(3200)
    return(__mm_ranks(x, w, method, mid, norm))
}

real colvector __mm_ranks(real colvector x, real colvector w,
 real scalar method, real scalar mid, real scalar norm) // sorted input assumed
{
    real scalar    i, I, i0, i1
    real colvector ranks, R
    pointer scalar wR
    pragma unset   R

    // compute running sum
    I = rows(x)
    if (I==0) return(J(0,1,.))
    if (rows(w)!=1) {
        // treat missing values as missing and use quad precision
        ranks = mm_colrunsum(w, 1, 1)
        if (ranks[I]>=.) return(J(I,1,.))
    }
    else ranks = (1::I) * w
    
    // no special treatment of ties; use fast computations 
    if (method==0 | method==4) {
        if (mid & norm) return((ranks :- w/2) / ranks[I])
        if (norm)       return(ranks / ranks[I])
        if (mid)        return(ranks :- w/2)
        return(ranks)
    }
    
    // normalize
    if (norm) ranks = ranks / ranks[I]
    
    // handle ties
    if (method==3) {            // use lowest rank
        for(i=I-1; i>=1; i--) {
            if (x[i]==x[i+1]) ranks[i] = ranks[i+1]
        }
    }
    else if (method==1) {       // use highest rank
        for(i=2; i<=I; i++) {
            if (x[i]==x[i-1]) ranks[i] = ranks[i-1]
        }
    }
    else if (method==2) {       // use average rank
        i0 = i1 = 1
        if (rows(w)==1) wR = &1
        else wR = &R
        for(i=2; i<=I; i++) {
            if (x[i0]==x[i]) {
                i1++
                if (i<I) continue
            }
            if (i0<i1) {
                R = (i0 \ i1)
                ranks[|R|] = J(i1-i0+1,1,1) :* mean(ranks[|R|], w[|*wR|])
            }
            i0 = i1 = i
        }
    }
    else _error(3498,"ties = " + strofreal(method) + " not allowed")
    
    // apply midpoint adjustment
    if (mid) {
        i0 = i1 = 0
        for (i=1; i<=I; i++) {
            if (ranks[i]!=i0) {
                i1 = (i0 + ranks[i])/2
                i0 = ranks[i]
            }
            ranks[i] = i1
        }
    }
    
    // done
    return(ranks)
}

end
