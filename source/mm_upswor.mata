*! version 1.0.1, Ben Jann, 12jun2007
version 9.1
mata:

// unequal probability sampling, without replacement
real colvector mm_upswor(real scalar n, real colvector w,
 | real scalar count, real scalar nowarn)
{
	real colvector p, ub, res
	real scalar i, j, nn, u, z

// check args
	if (args()<4) nowarn = 0
	if (args()<3) count = 0
	if (rows(w)<1) _error(3498, "no cases")
	if (n<0) _error(3300)
	if (n>=.) nn = rows(w)
	else nn = n
	if (nowarn==0) {
		if (rows(w)<nn) _error(3300, "n may not be larger than rows(w)")
	}

// no sampling
	if (rows(w)==1) {
		if (count) return(nn)
		return(J(nn,1,1))
	}

// draw sample
	p = mm_unorder2(rows(w)) //=> stable results even for very large N
	ub = w[p] / colsum(w) * nn
	if (nowarn==0) {
		z = colsum(ub:>1)
		if (z==1) _error(3300, "1 case has w_i*n/sum(w)>1")
		if (z>0) _error(3300, strofreal(z) + " cases have w_i*n/sum(w)>1")
	}
	ub = mm_colrunsum(ub)
	u = uniform(1,1)
	j=1
	if (count) {
		res = J(rows(w),1,0)
		for (i=1;i<=nn;i++) {
			while (u>ub[j]) j++
			res[j] = res[j]+1
			u = u + 1
		}
		return(res[invorder(p)])
	}
	res = J(nn,1,0)
	for (i=1;i<=nn;i++) {
		while (u>ub[j]) j++
		res[i] = p[j]
		u = u + 1
	}
	return(res)
}

end
