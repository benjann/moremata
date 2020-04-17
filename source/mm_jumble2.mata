*! version 1.0.0, Ben Jann, 29mar2006
version 9.0
mata:

transmorphic matrix mm_jumble2(transmorphic matrix x)
{
	return(x[mm_unorder2(rows(x)), .])
}

end
