*! version 1.0.4, Ben Jann, 22jun2006
version 9.2
mata:

real matrix mm_iqrange(real matrix X, | real colvector w,
 real scalar altdef)
{
	real matrix q

	if (args()<2) w = 1
	if (args()<3) altdef = 0
	q = mm_quantile(X, w, (.25 \ .75), altdef)
	return(q[2,]-q[1,])
}
end
