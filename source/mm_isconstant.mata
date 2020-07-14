*! version 1.0.1, Ben Jann, 13jul2020
version 9.0
mata:

real scalar mm_isconstant(X)
{
    if (length(X)<2) return(1)
    return(allof(X, X[1,1]))
}

end
