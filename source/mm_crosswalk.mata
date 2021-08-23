*! version 1.0.0  Ben Jann  23aug2021

version 11.2
mata:

transmorphic vector mm_crosswalk(
    transmorphic vector x,
    transmorphic vector from,
    transmorphic vector to,
    | transmorphic vector d, 
      transmorphic scalar n)
{
    real scalar         usehash, l, offset
    real rowvector      minmax
    transmorphic vector res

    // defaults
    if (args()<4) d = missingof(to)
    if (args()<5) n = 1e6
    
    // consistency checks
    if (eltype(x)!=eltype(from))  _error(3250, "{it:x} and {it:from} must be of same type")
    if (length(from)!=length(to)) _error(3200, "{it:from} and {it:to} must have same length")
    if (eltype(to)!=eltype(d))    _error(3250, "{it:to} and {it:d} must be of same type")
    if (length(d)!=1) {
        if (length(d)!=length(x)) _error(3200, "{it:x} and {it:d} not conformable")
    }
    
    // nothing to translate
    // - empty x
    if (!length(x)) return(J(rows(x), cols(x), missingof(to)))
    // - empty dictionary
    if (!length(to)) {
        if (length(d)!=1) return(J(1,1, rows(x)==rows(d) ? d : d'))
        return(J(rows(x), cols(x), d))
    }
    
    // pick algorithm
    usehash = 0
    if (!isreal(x)) usehash = 1
    else if (n<1)   usehash = 1
    else {
        if (n<.z) {
            if      (hasmissing(from))  usehash = 1
            else if (hasmissing(x))     usehash = 1
            else if (trunc(from)!=from) usehash = 1
            else if (trunc(x)!=x)       usehash = 1
        }
        if (!usehash) {
            minmax = minmax((minmax(from),minmax(x)))
            l = minmax[2] - minmax[1] + 1
            if (l>n) usehash = 1
        }
    }
    
    // use hash algorithm
    if (usehash) return(_mm_crosswalk_hash(x, from, to, d))
    
    // use index algorithm
    offset = minmax[1] - 1
    res = rows(to)!=1 ? J(l, 1, d[1]) : J(1, l, d[1])
    if (length(d)!=1) res[x :- offset] = (rows(res)!=1)==(rows(d)!=1) ? d : d'
    res[from :- offset] = to
    res = res[x :- offset]
    return(rows(res)==rows(x) ? res : res')
}

transmorphic vector mm_crosswalk_hash(
    transmorphic vector x,
    transmorphic vector from,
    transmorphic vector to,
    | transmorphic vector d)
{
    // defaults
    if (args()<4) d = missingof(to)
    
    // consistency checks
    if (eltype(x)!=eltype(from))  _error(3250, "{it:x} and {it:from} must be of same type")
    if (length(from)!=length(to)) _error(3200, "{it:from} and {it:to} must have same length")
    if (eltype(to)!=eltype(d))    _error(3250, "{it:to} and {it:d} must be of same type")
    if (length(d)!=1) {
        if (length(d)!=length(x)) _error(3200, "{it:x} and {it:d} not conformable")
    }
    
    // nothing to translate
    // - empty x
    if (!length(x)) return(J(rows(x), cols(x), missingof(to)))
    // - empty dictionary
    if (!length(to)) {
        if (length(d)!=1) return(J(1,1, rows(x)==rows(d) ? d : d'))
        return(J(rows(x), cols(x), d))
    }

    // run hash algorithm
    return(_mm_crosswalk_hash(x, from, to, d))
}

transmorphic vector _mm_crosswalk_hash(
    transmorphic vector x,
    transmorphic vector from,
    transmorphic vector to,
    transmorphic vector d)
{
    real scalar  i
    transmorphic A
    transmorphic scalar a
    real matrix notfound
    transmorphic vector res
    
    // build dictionary
    A = asarray_create(eltype(from))
    for (i=length(from); i; i--) asarray(A, from[i], to[i])
    
    // case 1: d is scalar
    if (length(d)==1) {
        res = J(rows(x), cols(x), missingof(to))
        asarray_notfound(A, d)
        for (i=length(x); i; i--) res[i] = asarray(A, x[i])
        return(res)
    }
    
    // case 2: d is non-scalar
    res = rows(x)==rows(d) ? d : d'
    notfound = J(0,0,.)
    asarray_notfound(A, notfound)
    for (i=length(x); i; i--) {
        a = asarray(A, x[i])
        if (a==notfound) continue // using asarray_contains() would be slower
        res[i] = a
    }
    return(res)
}

end
