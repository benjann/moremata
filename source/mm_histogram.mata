*! version 1.0.3, Ben Jann, 22jun2006
version 9.0
mata:

real matrix mm_histogram(
 real colvector x,     // data
 | real colvector w,  // weights
   real colvector g,  // grid (interval borders) (default: 10 intervals)
   real scalar dir)   // 0 right inclusive (default), else left inclusive
{
	real matrix H
	real scalar i

	if (args()<2) w = 1
	if (args()<3) g = rangen(min(x), max(x), 11)
	if (args()<4) dir = 0

	H = J(rows(g)-1, 3, 0)
	for (i=1; i<=rows(H); i++) {
		H[i,2] = g[i+1]-g[i]     // width of bins
		H[i,1] = g[i] + H[i,2]/2 // bin midpoints
	}
	H[.,3] = mm_exactbin(x, w, g, dir) :/ (mm_nobs(x,w) * H[.,2])
	return(H)
}
end
