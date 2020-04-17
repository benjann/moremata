*! version 1.0.0, Ben Jann, 07jun2006
version 9.2
mata:

real scalar mm_nobs(transmorphic matrix x, real colvector w)
{
	return(rows(w)==1 ? rows(x)*w :
	 (rows(x)==rows(w) ? colsum(w) : _error(3200)))
}
end
