*! version 1.1.0  09jul2020  Ben Jann
version 9.2
mata:

real scalar mm_nunique(transmorphic vector X)
{
    return(sum(mm_unique_tag(X)))
}

transmorphic vector mm_unique(transmorphic vector X, | real scalar order)
{
    if (rows(X)*cols(X)==0)  return(X)
    if (order==1 | order==2) return(select(X, mm_unique_tag(X, order)))
    // sorted
    if (cols(X)==1)    return(sort(select(X, mm_unique_tag(X)),1))
    if (!iscomplex(X)) return(sort(select(X, mm_unique_tag(X))',1)')
    return(transposeonly(sort(transposeonly(select(X, mm_unique_tag(X))),1)))
}

/*
real vector mm_unique_idx(transmorphic vector X, | real scalar order)
{
    if (rows(X)*cols(X)==0) return(J(rows(X),cols(X),.))
    if (cols(X)==1) return(select(1::rows(X), mm_unique_tag(X, order)))
                    return(select(1..cols(X), mm_unique_tag(X, order)))
}
*/

real vector mm_unique_tag(transmorphic vector X, | real scalar order)
{
    real vector p
    
    if (rows(X)*cols(X)==0) return(J(rows(X),cols(X),.))
    if (order==1) {
        // use stable sort order and tag first occurrence
        if (cols(X)==1)        p = mm_order(X,1,1)
        else if (iscomplex(X)) p = mm_order(transposeonly(X),1,1)'
        else                   p = mm_order(X',1,1)'
        p[p] = _mm_unique_tag(X[p])
    }
    else if (order==2) {
        // use stable sort order and tag last occurrence
        if (cols(X)==1)        p = mm_order(X,1,1)
        else if (iscomplex(X)) p = mm_order(transposeonly(X),1,1)'
        else                   p = mm_order(X',1,1)'
        p[p] = _mm_unique_tag(X[p], 1)
    }
    else {
        // use non-stable sort order (tag random element within ties)
        if (cols(X)==1)        p = order(X,1)
        else if (iscomplex(X)) p = order(transposeonly(X),1)'
        else                   p = order(X',1)'
        p[p] = _mm_unique_tag(X[p])
    }
    return(p)
}

real scalar _mm_nunique(transmorphic vector X)
{
    return(sum(_mm_unique_tag(X)))
}

transmorphic vector _mm_unique(transmorphic vector X)
{
    if (rows(X)*cols(X)==0) return(X)
    return(select(X, _mm_unique_tag(X)))
}

/*
real vector _mm_unique_idx(transmorphic vector X, | real scalar last)
{
    if (args()<2) last = 0
    if (rows(X)*cols(X)==0) return(J(rows(X),cols(X),.))
    if (cols(X)==1) return(select(1::rows(X), _mm_unique_tag(X, last)))
                    return(select(1..cols(X), _mm_unique_tag(X, last)))
}
*/

real vector _mm_unique_tag(transmorphic vector X, | real scalar last)
{
    real scalar r, c
    real vector p

    if (args()<2) last = 0
    r = rows(X); c = cols(X)
    if (r*c==0)      return(J(r,c,.))
    if (r==1 & c==1) return(J(1,1,1))
    if (last) {
        if (c==1) {
            p = (X :!= (X[|2\.|] \ X[1]))
            p[r] = 1 // should X be constant
            return(p)
        }
        p = (X :!= (X[|2\.|] , X[1]))
        p[c] = 1 // should X be constant
        return(p)
    }
    if (c==1) {
        p = (X :!= (X[r] \ X[|1\r-1|]))
        p[1] = 1 // should X be constant
        return(p)
    }
    p = (X :!= (X[c] , X[|1\c-1|]))
    p[1] = 1 // should X be constant
    return(p)
}

real scalar mm_nuniqrows(transmorphic matrix X)
{
    return(sum(mm_uniqrows_tag(X)))
}

transmorphic matrix mm_uniqrows(transmorphic matrix X, | real scalar order)
{
    if (rows(X)==0) return(X)
    if (order==1 | order==2) return(select(X, mm_uniqrows_tag(X, order)))
    return(mm_sort(select(X, mm_uniqrows_tag(X))))
}

/*
real vector mm_uniqrows_idx(transmorphic matrix X, | real scalar order)
{
    if (rows(X)==0) return(J(0,1,.))
    return(select(1::rows(X), mm_uniqrows_tag(X, order)))
}
*/

real vector mm_uniqrows_tag(transmorphic matrix X, | real scalar order)
{
    real vector p
    
    if (rows(X)==0) return(J(0,1,.))
    if (order==1) {
        p = mm_order(X,.,1)
        p[p] = _mm_uniqrows_tag(X[p,])
    }
    else if (order==2) {
        p = mm_order(X,.,1)
        p[p] = _mm_uniqrows_tag(X[p,], 1)
    }
    else {
        p = mm_order(X,.,0)
        p[p] = _mm_uniqrows_tag(X[p,])
    }
    return(p)
}

real scalar _mm_nuniqrows(transmorphic matrix X)
{
    return(sum(_mm_uniqrows_tag(X)))
}

transmorphic vector _mm_uniqrows(transmorphic matrix X)
{
    if (rows(X)==0) return(X)
    return(select(X, _mm_uniqrows_tag(X)))
}

/*
real vector _mm_uniqrows_idx(transmorphic matrix X, | real scalar last)
{
    if (args()<2) last = 0
    if (rows(X)==0) return(J(0,1,.))
    return(select(1::rows(X), _mm_uniqrows_tag(X, last)))
}
*/

real vector _mm_uniqrows_tag(transmorphic matrix X, | real scalar last)
{
    real scalar r, c
    real vector p

    if (args()<2) last = 0
    r = rows(X); c = cols(X)
    if (r==0) return(J(0,1,.))
    if (r==1) return(J(1,1,1))
    if (last) {
        if (c==0) p = J(r,1,0)
        else      p = (rowsum(X :!= (X[|2,1 \ .,.|] \ X[1,])):!=0)
        p[r] = 1 // should X be constant
        return(p)
    }
    if (c==0) p = J(r,1,0)
    else      p = (rowsum(X :!= (X[r,] \ X[|1,1 \ r-1,.|])):!=0)
    p[1] = 1 // should X be constant
    return(p)
}

end

