*! version 1.0.2, Ben Jann, 16may2007
version 9.1
mata:

real colvector mm_freq(transmorphic matrix x,
 | real colvector w, transmorphic matrix levels)
{
	real colvector p

	if (args()<2) w = 1
	if (args()<3) levels = .
	if (cols(x)==0) return(_mm_freq(x, w, levels))
	if (rows(w)==1) return(_mm_freq(sort(x,1..cols(x)), w, levels))
	p = order(x,1..cols(x))
	return(_mm_freq(x[p,], w[p,], levels))
}

real colvector _mm_freq(transmorphic matrix x,
 | real colvector w, transmorphic matrix levels)
{
	real scalar    i, j, l
	real colvector result

	if (args()<2) w = 1
	if (args()<3) levels = .
	if (rows(w)!=1 & rows(w)!=rows(x)) _error(3200)
	if (levels==.) levels = _mm_uniqrows(x)
	if (rows(x)==0) return(J(0,1, .))
	l = rows(levels)
	result = J(l,1,0)
	j = 1
	for (i=1; i<=rows(x); i++) {
		for (;j<=l;j++) {
			if (x[i,]==levels[j,]) break
		}
		if (j>l) break
		result[j] = result[j] + (rows(w)!=1 ? w[i] : w)
	}
	return(result)
}

real colvector mm_freq2(transmorphic matrix x,
 | real colvector w)
{
    real colvector p

    if (args()<2) w = 1
    if (cols(x)==0) return(_mm_freq2(x, w))
    p = order(x,1..cols(x))
    if (rows(w)==1) return(_mm_freq2(x[p,],w)[invorder(p)])
    return(_mm_freq2(x[p,],w[p,])[invorder(p)])
}

real colvector _mm_freq2(transmorphic matrix x,
 | real colvector w)
{
    real scalar    i, j
    real colvector result

    if (args()<2) w = 1
    if (rows(w)!=1 & rows(w)!=rows(x)) _error(3200)
    if (rows(x)==0) return(J(0, 1, .))
    result = J(rows(x),1,.)
    j = 1
    for (i=2; i<=rows(x); i++) {
        if (x[i,]!=x[i-1,]) {
            result[|j \ i-1|] = J(i-j, 1, (rows(w)==1 ?
                (i-j)*w : sum(w[|j \ i-1|])))
            j = i
        }
    }
    result[|j \ i-1|] = J(i-j, 1, (rows(w)==1 ?
        (i-j)*w : sum(w[|j \ i-1|])))
    return(result)
}

transmorphic matrix _mm_uniqrows( // uniqrows() for sorted X
 transmorphic matrix x)
{
        real scalar             i, j, n, ns
        transmorphic matrix     res

        if (rows(x)==0) return(J(0,cols(x), missingof(x)))
        if (cols(x)==0) return(J(1,0, missingof(x)))

        ns = 1
        n = rows(x)
        for (i=2;i<=n;i++) {
                if (x[i-1,]!=x[i,]) ns++
        }
        res = J(ns, cols(x), x[1,1])
        res[1,] = x[1,]
        for (i=j=2;i<=n;i++) {
                if (x[i-1,]!=x[i,]) res[j++,] = x[i,]
        }
        return(res)
}

end
