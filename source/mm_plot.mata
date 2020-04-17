*! version 1.0.0, Ben Jann, 03aug2006
version 9.2
mata:

void mm_plot(
 real matrix X,
 | string scalar type,
   string scalar opts)
{
	real scalar rc

	rc = _mm_plot(X, type, opts)
	if (rc) _error(rc)
}

real scalar _mm_plot(
 real matrix X,
 | string scalar type,
   string scalar opts)
{
	real scalar   n, N, k, j, rc, updata
	string scalar vars, plottype
	real vector   vid

	updata = st_updata()
	plottype = type
	if (plottype=="") plottype = "scatter"
	n = rows(X)
	k = cols(X)
	N = st_nobs()
	if (N<n & k>0) st_addobs(n-N)
	vid = J(1, k, .)
	for (j=1; j<=k; j++) {
		st_store((1,n), vid[j]=st_addvar("double", st_tempname()), X[,j])
		st_varlabel(vid[j], "Column " + strofreal(j))
		vars = vars + " " + st_varname(vid[j])
	}
	st_updata(updata)
	rc = _stata("twoway " + plottype + vars + ", " + opts)
	updata = st_updata() // twoway might change data-have-changed flag
	if (N<n & k>0) st_dropobsin((N+1,n))
	st_dropvar(vid)
	st_updata(updata)
	return(rc)
}

end
