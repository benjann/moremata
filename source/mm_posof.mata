*! version 1.0.0, Ben Jann, 24mar2006
version 9.0
mata:

real scalar mm_posof(transmorphic vector haystack, transmorphic scalar needle)
{
	real scalar i

	for (i=1; i<=length(haystack); i++) {
		if (haystack[i]==needle) return(i)
	}
	return(0)
}
end
