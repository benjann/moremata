*! version 2.0.0  08jul2020  Ben Jann
version 9.2
mata:

real colvector mm_ipolate(real colvector x0, real colvector y0,
    real colvector x1, | real scalar outer)
{
    real colvector p0, p1, y1
    
    if (args()<4) outer = 0
    if (rows(y0)!=rows(x0)) _error(3200)
    p0 = order(x0,1)
    p1 = order(x1,1)
    y1 = J(rows(x1),1,.)
    y1[p1] = _mm_ipolate(x0[p0], y0[p0], x1[p1], outer)
    return(y1)
}

real colvector _mm_ipolate(real colvector x0, real colvector y0,
    real colvector x1, | real scalar outer)
{   // x0 and x1 assumed sorted
    real scalar    n, i, ai, bi
    real colvector a, b, x, y

    if (args()<4) outer = 0
    n = rows(x0)
    if (rows(y0)!=n) _error(3200)
    if (n<=1) return(mm_fastipolate(x0, y0, x1, outer))
    // collapse data
    a = x0:!=(x0[n]\x0[|1\n-1|]) // tag first obs within panel
    a[1] = 1                     // (should x0 be constant)
    a = select(1::n, a)          // get index of first obs within panel
    if ((i=rows(a))==n) return(mm_fastipolate(x0, y0, x1, outer)) // no ties
    b = x0:!=(x0[|2\.|]\x0[1])   // tag last obs within panel
    b[n] = 1                     // (should x0 be constant)
    b = select(1::n, b)          // get index of last obs within panel
    x = y = J(i,1,.)
    for (;i;i--) {
        ai = a[i]; bi = b[i]
        x[i] = x0[ai]
        if (ai==bi) {
            y[i] = y0[ai]
            continue
        }
        y[i] = mean(y0[|ai\bi|],1)
    }
    return(mm_fastipolate(x, y, x1, outer))
}

real colvector mm_fastipolate(real colvector x0, real colvector y0,
    real colvector x1, | real scalar outer)
{   // x0 and x1 assumed sorted; x0 assumed unique
    real scalar    n, i, j
    real colvector y1

    if (args()<4) outer = 0
    j = rows(x0); n = i = rows(x1)
    if (rows(y0)!=j) _error(3200)
    if (j<1) return(J(n, 1, .))
    y1 = J(n,1,.)
    // handle values of x1 > max(x0)
    for (;i;i--) {
        if (x1[i]<=x0[j]) break
    }
    if (outer) {
        if (i<n) y1[|i+1\.|] = J(n-i, 1, y0[j])
    }
    // handle values of x1 within minmax(x0)
    for (;i;i--) {
        for (;j;j--) {
            if (x0[j]==x1[i]) {
                y1[i] = y0[j]
                break
            }
            if (x0[j]<x1[i]) {
                y1[i] = y0[j] + (y0[j+1]-y0[j]) * (x1[i]-x0[j])/(x0[j+1]-x0[j])
                break
            }
        }
        if (!j) break
    }
    // handle values of x1 < min(x0)
    if (outer) {
        if (i) y1[|1\i|] = J(i, 1, y0[1])
    }
    return(y1)
}

end
