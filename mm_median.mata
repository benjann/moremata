*! version 1.0.3, Ben Jann, 22jun2006
version 9.2
mata:

real rowvector mm_median(real matrix X, | real colvector w)
{
	if (args()==1) w = 1
	return(mm_quantile(X, w, .5))
}
end
