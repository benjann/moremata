*! version 1.0.0, Ben Jann, 01apr2006
version 9.1
mata:

// unequal probability sampling, with replacement
real colvector mm_upswr(real scalar n, real colvector w,
 | real scalar count)
{
	real colvector ub, res, u
	real scalar i, j, nn

// check args
	if (args()<3) count = 0
	if (rows(w)<1) _error(3498, "no cases")
	if (n>=.) nn = rows(w)
	else nn = n

// no sampling
	if (rows(w)==1) {
		if (count) return(nn)
		return(J(nn,1,1))
	}

// draw sample
	ub = mm_colrunsum(w)
	ub = ub/ub[rows(ub)]
	u = sort(uniform(nn,1),1)
	j=1
	if (count) {
		res = J(rows(w),1,0)
		for (i=1;i<=nn;i++) {
			while (u[i]>ub[j]) j++
			res[j] = res[j]+1
		}
		return(res)
	}
	res = J(nn,1,0)
	for (i=1;i<=nn;i++) {
		while (u[i]>ub[j]) j++
		res[i] = j
	}
	return(mm_jumble2(res)) //=> stable results even for very large n
}

end
