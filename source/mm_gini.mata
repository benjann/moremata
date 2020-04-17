*! version 1.0.3  26mar2008  Ben Jann
version 9.2
mata:

real rowvector mm_gini(real matrix X, | real colvector w)
{
    real rowvector result
    real scalar c, i

    if (args()==1) w = 1
    c = cols(X)
    if (c<1) return(J(1,0,.))
    if (c==1) return(_mm_gini(X, w))
    result = J(1, c, .)
    for (i=1; i<=c; i++) {
        result[1, i] = _mm_gini(X[,i], w)
    }
    return(result)
}

real scalar _mm_gini(real colvector x, real colvector w)
{
    real matrix mv

    if (rows(x)<1) return(.)
    mv = mm_meanvariance0((x, mm_ranks(x, w, 3, 1, 1)), w)
    return(mv[3,1] * 2 / mv[1,1])
}

end
