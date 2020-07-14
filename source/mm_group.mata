*! version 1.0.0  09jul2020  Ben Jann
version 9.2
mata:

real colvector mm_group(transmorphic matrix X, | real rowvector idx)
{
    real scalar    r, c
    real colvector p
    
    r = rows(X); c = cols(X)
    if (args()<2) {
        if (r<=1 | c==0) return(J(r,1,1))
        p = order(X, 1..c)
        p[p] = _mm_group(X[p,])
    }
    else {
        if (length(idx)>c) _error(3200)
        if (r<=1 | c==0) return(J(r,1,1))
        p = order(X, idx)
        p[p] = _mm_group(X[p, abs(idx)])
    }
    return(p)
}

real colvector _mm_group(transmorphic matrix X)
{
    real scalar r, c
    real colvector p
    
    r = rows(X); c = cols(X)
    if (r<=1 | c==0) return(J(r,1,1))
    if (c==1) p = (X :!= (X[r] \ X[|1\r-1|]))
    else      p = (rowsum(X :!= (X[r,] \ X[|1,1 \ r-1,.|])) :!= 0)
    return(mm_colrunsum(p))
}

end


