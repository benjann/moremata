*! version 1.0.0  29may2017  Christopher Baum and Ben Jann
version 9.2
mata:

real matrix mm_sqrt(real matrix A)
{
    real colvector  L
    real matrix     X
    pragma unset    L
    pragma unset    X
    
    if (isfleeting(A)) _symeigensystem(A, X, L)
    else                symeigensystem(A, X, L)
    return(X * diag(sqrt(L)) * X')
}

end
