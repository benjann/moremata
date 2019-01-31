*! version 1.0.2  16may2007
version 9.2
mata:

real scalar mm_nunique(transmorphic vector x)
{
	if (cols(x)==1) return(mm_npanels(sort(x,1)))
	else {
		if (iscomplex(x))
		  return(mm_npanels(sort(transposeonly(x),1)))
		else
		  return(mm_npanels(sort(x',1)))
	}
}

real scalar mm_nuniqrows(transmorphic matrix x)
{
	real scalar i, n, ns
	real colvector p

	if (rows(x)==0) return(0)
	if (cols(x)==0) return(1)

	p = order(x, 1..cols(x))
	ns = 1
	n = rows(x)
	for (i=2;i<=n;i++) {
		if (x[p[i-1],]!=x[p[i],]) ns++
	}
	return(ns)
}

end
