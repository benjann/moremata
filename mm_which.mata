*! version 1.0.2, Ben Jann, 17apr2007
version 9.2
mata:

real matrix mm_which(real vector I)
{
	if (cols(I)!=1) return(select(1..cols(I), I))
	else return(select(1::rows(I), I))
}

end
