*! version 1.0.1, Ben Jann, 22jun2006
version 9.0
mata:

real rowvector mm_colvar(real matrix X, |real colvector w)
{
	real rowvector result
	real scalar j, k

	if (args()==1) w = 1

	result = J(1, k=cols(X), .)
	for (j=1; j<=k; j++) {
		result[j] = variance(X[,j], w)
	}
	return(result)
}

real matrix mm_meancolvar(real matrix X, |real colvector w)
{
	real matrix result
	real scalar j, k

	if (args()==1) w = 1

	result = J(2, k=cols(X), .)
	for (j=1; j<=k; j++) {
		result[.,j] = meanvariance(X[,j], w)
	}
	return(result)
}

end
