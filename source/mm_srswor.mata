*! version 1.0.0, Ben Jann, 29mar2006
version 9.1
mata:

// simple random sample, without replacement
real colvector mm_srswor(real scalar n, real scalar N,
 | real scalar count)
{
	real colvector res, u
	real scalar i, nn

// check args
	if (args()<3) count = 0
	if (N>=.) _error(3351)
	if (N<1) _error(3300)
	if (n<0) _error(3300)
	if (n>=.) nn = N
	else nn = n
	if (N<nn) _error(3300, "n may not be larger than N")

// no sampling
	if (N==1) {
		if (count) return(nn)
		return(J(nn,1,1))
	}
	if (nn==0) {
		if (count) return(J(N,1,0))
		return(J(0,1,.))
	}

// draw sample
	u = mm_unorder2(N)[|1 \ nn|] //=> stable results even for very large N
	if (count==0) return(u)
	res = J(N,1,0)
	for (i=1;i<=nn;i++) {
		res[u[i]] = 1
	}
	return(res)
}
end
