*! version 1.0.0, Ben Jann, 16mar2006
version 9.1
mata:

real colvector mm_cut(real colvector x, real vector at, | real scalar sorted)
{
	real scalar i, j
	real colvector result, p

	result = J(rows(x),1,.)
	j = length(at)
	if (args()<3) sorted = 0
	if (sorted) {
		for (i=rows(x); i>0; i--) {
			if (x[i]>=.) continue
			for (; j>0; j--) {
				if (at[j]<=x[i]) break
			}
			if (j>0) result[i,1] = at[j]
		}
		return(result)
	}
	p = order(x,1)
	for (i=rows(x); i>0; i--) {
		if (x[p[i]]>=.) continue
		for (; j>0; j--) {
			if (at[j]<=x[p[i]]) break
		}
		if (j>0) result[p[i],1] = at[j]
	}
	return(result)
}
end
