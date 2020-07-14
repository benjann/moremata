*! version 1.1.0  10jul2020  Ben Jann
version 9.2
mata:

real matrix mm_ecdf(real matrix X, | real colvector w, real scalar mid, 
    real scalar nonorm, real scalar brk)
{
    if (args()<2) w = 1
    if (args()<3) mid = 0
    if (args()<4) nonorm = 0
    if (args()<5) brk = 0
    return(mm_ranks(X, w, (brk!=0 ? 0 : 3), mid, (nonorm!=0 ? 0 : 1)))
}

real colvector _mm_ecdf(real colvector X, | real colvector w, real scalar mid, 
    real scalar nonorm, real scalar brk)
{
    if (args()<2) w = 1
    if (args()<3) mid = 0
    if (args()<4) nonorm = 0
    if (args()<5) brk = 0
    return(_mm_ranks(X, w, (brk!=0 ? 0 : 3), mid, (nonorm!=0 ? 0 : 1)))
}

real matrix mm_ecdf2(real colvector X, | real colvector w, real scalar mid, 
    real scalar nonorm)
{
    real colvector p

    if (args()<2) w = 1
    if (args()<3) mid = 0
    if (args()<4) nonorm = 0
    if (rows(X)==0) return(J(0,2,.))
    if (rows(w)==1) {
        p = order(X,1)
        return(_mm_ecdf2(X[p], w, mid, nonorm))
    }
    if (rows(w)!=rows(X)) _error(3200)
    p = order((X,w),(1,2))
    return(_mm_ecdf2(X[p], w[p], mid, nonorm))
}

real matrix _mm_ecdf2(real colvector X, | real colvector w, real scalar mid, 
    real scalar nonorm)
{
    real scalar    n
    real colvector p, cdf

    if (args()<2) w = 1
    if (args()<3) mid = 0
    if (args()<4) nonorm = 0

    // compute running sum
    n = rows(X)
    if (rows(w)==1) cdf = (1::n) * w
    else {
        if (rows(w)!=n) _error(3200)
        // treat missing values as missing and use quad precision
        cdf = mm_colrunsum(w, 1, 1)
        if (cdf[n]>=.) cdf = J(n,1,.)
    }

    // remove ties (select last obs in each group)
    p = select(1::n, _mm_unique_tag(X, 1))
    cdf = cdf[p]
    
    // normalize
    if (!nonorm) cdf = cdf / cdf[rows(cdf)]
    
    // midrank adjustment
    if (mid) cdf = cdf :- mm_diff(0\cdf)/2
    
    // return result
    return(X[p], cdf)
}

end
