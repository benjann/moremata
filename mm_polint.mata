*! version 1.0.0, Ben Jann, 12jul2006
version 9.2
mata:

real vector mm_polint(
 real vector xa,
 real vector ya,
 real vector x,
 | real scalar degree) // degree
{
	real scalar  n, m, i, j, a, b
	real vector  y

	if ( args()<4 ) degree = 1
	n = length(xa)
	if (length(ya)!=n) _error(3200)
	if ( degree<1 ) _error(3498, "degree must be 1 or larger")
	if ( degree>=n ) _error(3498, "degree must be smaller than length of data")
	if ( degree!=trunc(degree) ) _error(3498, "degree must be integer")

	m = degree + 1
	y = J(rows(x), cols(x), .)
	j = 0
	for (i=1; i<=length(x); i++) {
		mm_hunt(xa, x[i], j)
		a = min((max((trunc(j-(m-1)/2+.5), 1)), n+1-m))
		b = a+m-1
		y[i] = _mm_polint(xa[|a \ b|], ya[|a \ b|], x[i])
	}
	return(y)
}

real scalar _mm_polint(
// translation of -polint- from Press et al. (1992), Numerical
// Recipes in C, p. 109-110, http://www.numerical-recipes.com
 real vector xa,
 real vector ya,
 real scalar x)
{
	real scalar  i, m, n, ns
	real scalar  y, den, dif, dift, ho, hp, w
	real vector  c, d

	n = length(xa)
	if ( length(ya) != n ) _error(3200)

	ns = 1
	dif = abs(x-xa[1])
	c = J(1, n, .)
	d = J(1, n, .)
	for (i=1; i<=n; i++) {
		if ( (dift=abs(x-xa[i])) < dif) {
			ns  = i
			dif = dift
		}
		c[i] = ya[i]
		d[i] = ya[i]
	}
	y = ya[ns--]
	for (m=1; m<n; m++) {
		for (i=1; i<=n-m; i++) {
			ho = xa[i] - x
			hp = xa[i+m] - x
			w = c[i+1] - d[i]
			if ( (den=ho-hp) == 0.0)
			 _error(3498, "invalid input data") // This error can occur
			 // only if two input xa's are (to within roundoff) identical.
			den = w / den
			d[i] = hp * den
			c[i] = ho * den
		}
		if ( 2*ns < (n-m) ) y = y + c[ns+1]
		else y = y + d[ns--]
	}
	return(y)
}

end
