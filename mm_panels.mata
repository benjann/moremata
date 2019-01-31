*! version 1.0.3, Ben Jann, 03may2007
version 9.2
mata:

void function mm_panels(
 transmorphic vector X1,       // upper level ID variable (e.g. strata)
 transmorphic matrix info1,    // will be replaced (unless X1 and X2 missing)
 | transmorphic vector X2,     // lower level ID variable (e.g. clusters)
   transmorphic matrix info2)  // will be replaced (unless X2 missing)
{
	real scalar i, j, nc, b, e, b2, e2

	if (length(X1)>0 & X1!=.) info1 = _mm_panels(X1)
	if (args()<3) return
	if (length(X2)==0 | X2==.) {
		if (length(X1)>0 & X1!=.) info1 = info1, info1
		return
	}
	if (length(X1)==0 | X1==.) {
		info2 = _mm_panels(X2)
		info1 = colsum(info2), rows(info2)
		return
	}
	if (length(X1)!=length(X2)) _error(3200)
	info1 = info1, J(rows(info1),1,.)
	e = 0
	for (i=1; i<=rows(info1); i++) {
		nc = 1
		b = e + 2
		e = e + info1[i,1]
		for (j=b; j<=e; j++) {
			if (X2[j]!=X2[j-1]) nc++
		}
		info1[i,2] = nc
	}
	if (args()<4) return
	info2 = J(colsum(info1[.,2]), 1, .)
	e = e2 = 0
	for (i=1; i<=rows(info1); i++) {
		b = e + 1
		e = e + info1[i,1]
		b2 = e2 + 1
		e2 = e2 + info1[i,2]
		info2[|b2 \ e2|] = _mm_panels(X2[|b \ e|], info1[i,2])
	}
}

real colvector _mm_panels(transmorphic vector X, | real scalar np)
{
	real scalar i, j, r
	real colvector res

	r = length(X)
	if (r<1) return(J(0,1,.))
	if (args()<2) np = r
	res = J(np, 1, 1)
	j = 1
	for (i=2; i<=r; i++) {
		if (X[i]==X[i-1]) res[j] = res[j] + 1
		else j++
	}
	if (j==r) return(res)
	return(res[|1 \ j|])
}

// /*Old (slow) version of _mm_panels()*/
//real colvector _mm_panels(transmorphic vector X, | real scalar np)
//{
//	real scalar i, j, n
//	real colvector res
//
//	if (args()<2) np = mm_npanels(X)
//	if (length(X)==0) return(J(0,1,.))
//	res = J(np, 1, .)
//	n = j = 1
//	for (i=2; i<=length(X); i++) {
//		if (X[i]!=X[i-1]) {
//			res[j++] = n
//			n = 1
//		}
//		else n++
//	}
//	res[j] = n
//	return(res)
//}

real scalar mm_npanels(vector X)
{
	real scalar i, np

	if (length(X)==0) return(0)
	np = 1
	for (i=2; i<=length(X); i++) {
		if (X[i]!=X[i-1]) np++
	}
	return(np)
}

end
