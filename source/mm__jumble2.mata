*! version 1.0.0, Ben Jann, 29mar2006
version 9.0
mata:

function mm__jumble2(transmorphic matrix x)
{
	_collate(x,mm_unorder2(rows(x)))
}

end
