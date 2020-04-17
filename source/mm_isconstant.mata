*! version 1.0.0, Ben Jann, 01mar2007
version 9.0
mata:

real scalar mm_isconstant(X)
{
	if (length(X)<2) return(1)
	return(all(X:==X[1,1]))
}

end
