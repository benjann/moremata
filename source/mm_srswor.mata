*! version 1.0.1, Ben Jann, 23oct2020
version 9.1
mata:

// simple random sample, without replacement
real colvector mm_srswor(real scalar n, real scalar N,
 | real scalar count, real scalar alt)
{
	real colvector res, u
	real scalar i, nn

// check args
	if (args()<3) count = 0
	if (args()<4) alt = 0
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
	if (alt) {
		if (nn^2 < N/2) u = _mm_srswor_a(n, N)
		else            u = _mm_srswor_b(n, N)
	}
	else u = mm_unorder2(N)[|1 \ nn|] //=> stable results even for very large N
	if (count==0) return(u)
	res = J(N,1,0)
	for (i=1;i<=nn;i++) {
		res[u[i]] = 1
	}
	return(res)
}

real colvector _mm_srswor_a(real scalar n, real scalar N)
{   // naive algorithm; fast when n is small compared to N
    real scalar    i, u
    real colvector H
    
    H = J(n, 1, .)
    H[1] = ceil(uniform(1,1)*N)
    for (i=2; i<=n; i++) {
        while (anyof(H[|1 \ i|], u = ceil(uniform(1,1)*N))) continue
        H[i] = u
    }
    return(H)
}

real colvector _mm_srswor_b(real scalar n, real scalar N)
{   // Fisher–Yates shuffle; fast when n is large compared to N
    real scalar    i, j, k
    real colvector I, H
    
    k = N
    I = 1::k
    H = J(n, 1, .)
    for (i=1; i<=n; i++) {
        j = ceil(uniform(1,1)*k)
        H[i] = I[j]
        I[j] = I[k]
        k--
    }
    return(H)
}


end
