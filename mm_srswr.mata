*! version 1.0.0, Ben Jann, 27mar2006
version 9.1
mata:

// simple random sample, with replacement
real colvector mm_srswr(real scalar n, real scalar N,
 | real scalar count)
{
	real colvector res, u
	real scalar i, nn

// check args
	if (args()<3) count = 0
	if (N>=.) _error(3351)
	if (N<1) _error(3300)
	if (n>=.) nn = N
	else nn = n

// no sampling
	if (N==1) {
		if (count) return(nn)
		return(J(nn,1,1))
	}

// draw sample
	u = ceil(uniform(nn,1)*N)
	if (count==0) return(u)
	res = J(N,1,0)
	for (i=1;i<=nn;i++) {
		res[u[i]] = res[u[i]] + 1
	}
	return(res)
}

end
