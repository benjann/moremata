*! version 1.0.2, Ben Jann, 24mar2006
version 9.0
mata:

real matrix mm_realofstr(string matrix S)
{
	real matrix R
	string scalar tmp
	real scalar i, j

	R = J(rows(S), cols(S), .)
	tmp = st_tempname()
	for (i=1; i<=rows(S); i++) {
		for (j=1; j<=cols(S); j++) {
			stata("scalar " + tmp + "=real(`" + `"""' + S[i,j] + `"""' + "')")
			R[i,j] = st_numscalar(tmp)
		}
	}
	stata("capture scalar drop " + tmp)
	return(R)
}
end
