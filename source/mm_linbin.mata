*! version 1.0.9  12aug2020  Ben Jann
version 9.0
mata:

real colvector mm_linbin(
    real colvector x,  // data points
    real colvector w,  // weights
    real colvector g)  // grid points (assumed sorted)
{
    real colvector  p

    p = mm_order(x, ., 1) // stable sort order
    if (rows(w)==1) return(__mm_linbin(x[p], g)*w)
                    return(__mm_linbin_w(x[p], w[p], g))
}

real colvector _mm_linbin(
    real colvector x,  // data points (assumed sorted)
    real colvector w,  // weights
    real colvector g)  // grid points (assumed sorted)
{
    if (rows(w)==1) return(__mm_linbin(x, g)*w)
                    return(__mm_linbin_w(x, w, g))
}

real colvector __mm_linbin(
    real colvector x,  // data points (assumed sorted)
    real colvector g)  // grid points (assumed sorted)
{
    real scalar     i, j, g1, g0, i1, i0
    real colvector  c, xx

    j = rows(g)
    i = rows(x)
    if (i<1) return(J(j, 1, 0))
    if (j<2) return(i)
    c = J(j, 1, 0)
    // data above grid
    g0 = g[j]
    i1 = i
    for (;i;i--) {
        if (x[i]<g0) break
    }
    i0 = i+1
    c[j] = (i1-i0+1)
    // data within grid range
    g1 = g0
    for (j--; j; j--) {
        g1 = g0
        g0 = g[j]
        i1 = i
        for (;i;i--) {
            if (x[i]<g0) break
        }
        i0 = i+1
        if (i1<i0) continue // no data in current bin
        xx = x[|i0\i1|]
        c[j] =            sum(g1 :- xx) / (g1 - g0)
        c[j+1] = c[j+1] + sum(xx :- g0) / (g1 - g0)
    }
    // data below grid
    if (i) c[1] = c[1] + i
    return(c)
}

real colvector __mm_linbin_w(
    real colvector x,  // data points (assumed sorted)
    real colvector w,  // weights
    real colvector g)  // grid points (assumed sorted)
{
    real scalar     i, j, g1, g0, i1, i0
    real colvector  c, xx, ww

    j = rows(g)
    i = rows(x)
    if (rows(w)!=i) _error(3200)
    if (i<1) return(J(j, 1, 0))
    if (j<2) return(sum(w))
    c = J(j, 1, 0)
    // data above grid
    g0 = g[j]
    i1 = i
    for (;i;i--) {
        if (x[i]<g0) break
    }
    i0 = i+1
    if (i0<=i1) c[j] = sum(w[|i0 \ i1|])
    // data within grid range
    g1 = g0
    for (j--; j; j--) {
        g1 = g0
        g0 = g[j]
        i1 = i
        for (;i;i--) {
            if (x[i]<g0) break
        }
        i0 = i+1
        if (i1<i0) continue // no data in current bin
        xx = x[|i0\i1|]
        ww = w[|i0\i1|]
        c[j] =            sum(ww:*(g1 :- xx)) / (g1 - g0)
        c[j+1] = c[j+1] + sum(ww:*(xx :- g0)) / (g1 - g0)
    }
    // data below grid
    if (i) c[1] = c[1] + sum(w[|1 \ i|])
    return(c)
}

real colvector mm_fastlinbin(
    real colvector x,  // data points
    real colvector w,  // weights
    real colvector g)  // grid points (assumed sorted)
{
    real scalar     g1, i, j, d, z, ng
    real colvector  c

    ng = rows(g)
    i  = rows(x)
    if (ng<1) _error(3200)
    if (i<1)  return(J(ng, 1, 0))
    if (ng<2) return(mm_nobs(x,w))
    if (rows(w)!=1) return(_mm_fastlinbin_w(x, w, g))
    g1 = g[1]
    d = (g[ng] - g1)/(ng-1)
    c = J(ng, 1, 0)
    for (; i; i--) {
        z = (x[i] - g1) / d
        if (z<=0) {
            c[1] = c[1] +  1
            continue
        }
        if (z>=(ng-1)) {
            c[ng] = c[ng] +  1
            continue
        }
        j = trunc(z) + 1
        c[j] = c[j] + (j-z)
        j++
        c[j] = c[j] + (2-j+z)
    }
    return(c*w)
}

real colvector _mm_fastlinbin_w(
    real colvector x,  // data points
    real colvector w,  // weights
    real colvector g)  // grid points (assumed sorted)
{
    real scalar     g1, i, j, d, z, ng
    real colvector  c

    ng = rows(g)
    i  = rows(x)
    if (rows(w)!=i) _error(3200)
    g1 = g[1]
    d = (g[ng] - g1)/(ng-1)
    c = J(ng, 1, 0)
    for (; i; i--) {
        z = (x[i] - g1) / d
        if (z<=0) {
            c[1] = c[1] +  w[i]
            continue
        }
        if (z>=(ng-1)) {
            c[ng] = c[ng] +  w[i]
            continue
        }
        j = trunc(z) + 1
        c[j] = c[j] + w[i] * (j-z)
        j++
        c[j] = c[j] + w[i] * (2-j+z)
    }
    return(c)
}

end
