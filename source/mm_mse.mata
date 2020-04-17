*! version 1.0.0, Ben Jann, 22jun2006
version 9.0
mata:

real matrix mm_mse(real matrix X, real colvector w, real rowvector mu)
{
	real rowvector CP

	CP = cross(w,0, X,1)
	if (missing(CP) | CP[cols(CP)]==0) return(J(cols(X), cols(X), .))
	return(mm_sse(X, w, mu) :/ CP[cols(CP)])
}

real matrix mm_sse(real matrix X, real colvector w, real rowvector mu)
{
	real rowvector CP

	CP = cross(w,0, X,1)
	if (missing(CP) | CP[cols(CP)]==0) return(J(cols(X), cols(X), .))
	return(crossdev(X,0,mu, w, X,0,mu))
}

real rowvector mm_colmse(real matrix X, real colvector w, real rowvector mu)
{
	real rowvector result
	real scalar j, k

	result = J(1, k=cols(X), .)
	for (j=1; j<=k; j++) {
		result[j] = mm_mse(X[,j], w, mu[j])
	}
	return(result)
}

real rowvector mm_colsse(real matrix X, real colvector w, real rowvector mu)
{
	real rowvector result
	real scalar j, k

	result = J(1, k=cols(X), .)
	for (j=1; j<=k; j++) {
		result[j] = mm_sse(X[,j], w, mu[j])
	}
	return(result)
}

end
