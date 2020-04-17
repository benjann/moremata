*! version 1.0.4, Ben Jann, 26jun2006
version 9.0
mata:

real colvector mm_makegrid(
 real colvector x, // data points
 | real scalar M,  // number of bins (default: 401)
   real scalar e,  // extend range by e on each side
   real scalar min,  // min
   real scalar max)  // max
{
	real scalar a, b

	a = (min<. ? min : min(x) - (e<. ? e : 0))
	b = (max<. ? max : max(x) + (e<. ? e : 0))
	return( rangen(a, b, (M<. ? M : 512)) )
}

end
