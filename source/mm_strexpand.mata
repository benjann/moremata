*! version 1.0.4, Ben Jann, 30apr2007
version 9.0
mata:

string scalar mm_strexpand(string scalar s, string vector slist,
 | string scalar def, real scalar unique, string scalar errtxt)
{
	real scalar   err
	string scalar res

	if (args()<5) errtxt = `"""' + s + `"" invalid"'
	if (args()<4) unique = 0
	err = _mm_strexpand(res, s, slist, def, unique)
	if (err) _error(err, errtxt)
	return(res)
}

real scalar _mm_strexpand(res, string scalar s,
 string vector slist, | string scalar def, real scalar unique)
{
	real scalar i, l, match

	if (s=="") {
		res = def
		return(0)
	}
	if (args()<5) unique = 0
	l = strlen(s)
	if (unique) {
		match = 0
		for (i=1; i<=length(slist); i++) {
			if (s==substr(slist[i], 1, l)) {
				if (match) return(3498)
				match = i
			}
		}
		if (match) {
			res = slist[match]
			return(0)
		}
	}
	else {
		for (i=1; i<=length(slist); i++) {
			if (s==substr(slist[i], 1, l)) {
				res = slist[i]
				return(0)
			}
		}
	}
	return(3499)
}
end
