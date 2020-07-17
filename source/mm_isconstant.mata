*! version 1.0.2, Ben Jann, 16jul2020
version 9.0
mata:

real scalar mm_isconstant(transmorphic matrix X)
{
    if (length(X)<2) return(1)
    return(allof(X, X[1,1]))
}

real scalar mm_issorted(transmorphic vector x, | real scalar descending)
{
    if (ispointer(x)) _error(3250)
    if (length(x)<=1) return(1)
    if (args()<2) descending = 0
    if (descending) {
        if (isstring(x)) {
            if (cols(x)==1) return(all(x[|1 \ length(x)-1|] :>= x[|2 \ .|]))
                            return(all(x[|1 \ length(x)-1|] :>= x[|2 \ .|]))
        }
        if (cols(x)==1) return(all((.z \ x[|1 \ length(x)-1|]) :>= x))
                        return(all((.z, x[|1 \ length(x)-1|]) :>= x))
    }
    if (isstring(x)) {
        if (cols(x)==1) return(all(x[|1 \ length(x)-1|] :<= x[|2 \ .|]))
                        return(all(x[|1 \ length(x)-1|] :<= x[|2 \ .|]))
    }
    if (cols(x)==1) return(all(x :<= (x[|2 \ .|] \ .z)))
                    return(all(x :<= (x[|2 \ .|], .z)))
}

end
