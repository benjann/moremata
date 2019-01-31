*! version 1.0.6  03aug2007  Ben Jann
version 9.2
mata:

real matrix mm_ranks(real matrix X, | real colvector w,
 real scalar method, real scalar mid, real scalar norm)
{
    real rowvector result
    real scalar c, r, i

    if (args()<2) w = 1
    if (args()<3) method = 0
    if (args()<4) mid = 0
    if (args()<5) norm = 0
    c = cols(X)
    r = rows(X)
    if (rows(w)!=1 & rows(w)!=r) _error(3200)
    if (c<1) return(J(r,0,.))
    if (c==1) return(_mm_ranks(X, w, method, mid, norm))
    result = J(r, c, .)
    for (i=1; i<=c; i++) {
        result[,i] = _mm_ranks(X[,i], w, method, mid, norm)
    }
    return(result)
}

real colvector _mm_ranks(real colvector x, real colvector w,
 real scalar method, real scalar mid, real scalar norm)
{
    real scalar    I
    real colvector p

    I = rows(x)
    if (I==0) return(J(0,1,.))
    if (method==4 & rows(w)!=1) p = order((x,w),(1,2)) // -set sortseed- for reproducibility
    else if (method==0 | method==4) p = order(x,1) // -set sortseed- for reproducibility
    else if (method==1 & rows(w)!=1) p = order((x,w,(1::I)),(1,2,3)) // stable sort order
    else p = order((x,(1::I)),(1,2)) // stable sort order
    return(__mm_ranks(x[p],(rows(w)!=1 ? w[p] : w), method, mid, norm)[invorder(p)])
}

real colvector __mm_ranks(real colvector x, real colvector w,
 real scalar method, real scalar mid, real scalar norm)   // sorted input assumed
{
    return(_mm_ranks_mid(x, ___mm_ranks(x, w, method, norm), mid))
}

real colvector ___mm_ranks(real colvector x, real colvector w,
 real scalar method, real scalar norm)                    // sorted input assumed
{
    real scalar    i, I, i0, i1
    real colvector ranks, R
    pointer scalar wR

    I = rows(x)
    if (I==0) return(J(0,1,.))
    if (rows(w)!=1) {
        if (missing(w)) return(J(I,1,.))
        ranks = mm_colrunsum(w)
    }
    else ranks = (1::I) * w
    if (norm) ranks = ranks / ranks[I]
    if (method==0 | method==4) return(ranks)    // no treatment of ties
    if (method==3) {                            // ties lowest
        for(i=I-1; i>=1; i--) {
            if (x[i]==x[i+1]) ranks[i] = ranks[i+1]
        }
        if (method==3) return(ranks)
    }
    if (method==1) {                            // ties highest
        for(i=2; i<=I; i++) {
            if (x[i]==x[i-1]) ranks[i] = ranks[i-1]
        }
        return(ranks)
    }
    if (method!=2) _error(3498,"ties = " + strofreal(method) + " not allowed")
    i0 = i1 = 1                                 // ties mean
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
    return(ranks)
}

real colvector _mm_ranks_mid(real colvector x, real colvector ranks,
    real scalar mid)
{
    real scalar     r, lr, i

    if (rows(ranks)<1 | mid==0) return(ranks)
    r = ranks[1]
    ranks[1] =  r/2
    lr = r
    for(i=2; i<=rows(ranks); i++) {
        if (x[i]==x[i-1]) ranks[i] = ranks[i-1]
        else {
            r = ranks[i]
            ranks[i] = lr + (r-lr)/2
            lr = r
        }
    }
    return(ranks)
}
end
