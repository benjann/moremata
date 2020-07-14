*! version 1.0.0  09jul2020  Ben Jann
version 9.2
mata:

transmorphic matrix mm_sort(transmorphic matrix X,
    | real rowvector idx, real scalar stable)
{
    if (args()<2) return(X[mm_order(X), .])
    if (args()<3) return(X[mm_order(X, idx), .])
    return(X[mm_order(X, idx, stable), .])
}

real colvector mm_order(transmorphic matrix X,
    | real rowvector idx, real scalar stable)
{
    real scalar r, c, l
    
    if (args()<3) stable = 0
    // redirect to order()
    c = cols(X)
    if (!stable) {
        if (args()>=2 & idx!=.) return(order(X, idx))
        return(order(X, (c>0 ? 1..c : J(1,0,.))))
    }
    // check idx
    l = (idx==. ? 0 : length(idx))
    if (l>c) _error(3200)
    // one row or less: nothing to do
    r = rows(X)
    if (r<=1) return(J(r,1,1))
    // zero columns: return index
    if (c==0) return(1::r)
    // X is real: append index for stable order
    if (isreal(X)) {
        if (l==0) return(order((X,(1::r)), 1..c+1))
        if (anyof(idx,0)) _error(3300)
        return(order((X[,abs(idx)],(1::r)), ((1..l):*sign(idx), l+1)))
    }
    // X is not real; need to create a group id
    if (l==0) return(order((mm_group(X), (1::r)), (1,2)))
    return(order((mm_group(X, idx), (1::r)), (1,2)))
}

end

