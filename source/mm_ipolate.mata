*! version 1.0.6, Ben Jann, 10jul2006
version 9.2
mata:

real colvector mm_ipolate(real colvector x, real colvector y,
 real colvector xnew, | real scalar outer)
{
	real scalar i, j0b, j0e, j1b, j1e, r, xlo, xup, xi, y0, y1
	real colvector p, pnew, ynew

	r = rows(x)
	if (rows(y)!=r) _error(3200)
	if (r<1) return(J(rows(xnew), 1, .))
	if (args()<4) outer = 0

	p = order(x, 1)
	pnew = order(xnew, 1)
	ynew = J(rows(xnew),1,.)
	xlo = x[p[1]]; xup = x[p[rows(p)]]
	j0b = j1e = j0e = 1
	for (i=1; i<=rows(xnew); i++) {
		xi = xnew[pnew[i]]
		if (outer==0) {
			if (xi<xlo) continue
			if (xi>xup) return(ynew)
		}
		while (j0e<r) {
			if (x[p[j0e+1]]>xi) break
			j0e++
			if (x[p[j0e]]>x[p[j0b]]) j0b = j0e
		}
		if (j0e>=j1e) {
			j1b = j0e
			while (j1b<=r) {
				if (x[p[j1b]]>=xi) break
				j1b++
			}
			if (j1b>r) j1b = r
			j1e = j1b
			while (j1e<r) {
				if (x[p[j1e+1]]>x[p[j1b]]) break
				j1e++
			}
		}
		y0 = (j0b==j0e ? y[p[j0b]] : mean(y[p[|j0b \ j0e|]],1))
		y1 = (j1b==j1e ? y[p[j1b]] : mean(y[p[|j1b \ j1e|]],1))
		if (outer) {
			if (xi<xlo) {
				ynew[pnew[i]] = y1
				continue
			}
			if (xi>xup) {
				ynew[pnew[|i \ rows(pnew)|]] = J(rows(pnew)-i+1,1,y0)
				return(ynew)
			}
		}
		ynew[pnew[i]] = ( j0e==j1e ? y0 :
		 y0 + (y1-y0) *  (xi-x[p[j0e]])/(x[p[j1b]]-x[p[j0e]]) )
	}
	return(ynew)
}
end
