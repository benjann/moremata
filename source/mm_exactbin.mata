*! version 1.0.6  12aug2020  Ben Jann
version 9.0
mata:

real colvector mm_exactbin(
    real colvector x,   // data
    real colvector w,   // weights
    real colvector g,   // grid (assumed sorted)
    | real scalar dir,  // 0 left inclusive (default),  else right inclusive
      real scalar include) // include data outside grid in first/last bin
{
    real colvector  p

    if (args()<4) dir = 0
    if (args()<5) include = 0
    p = mm_order(x, ., 1) // stable sort order
    if (rows(w)==1) return(__mm_exactbin(x[p], g, dir, include)*w)
                    return(__mm_exactbin_w(x[p], w[p], g, dir, include))
}

real colvector _mm_exactbin(
    real colvector x,   // data (assumed sorted)
    real colvector w,   // weights
    real colvector g,   // grid (assumed sorted)
    | real scalar dir,  // 0 left inclusive (default),  else right inclusive
      real scalar include) // include data outside grid in first/last bin
{
    if (args()<4) dir = 0
    if (args()<5) include = 0
    if (rows(w)==1) return(__mm_exactbin(x, g, dir, include)*w)
                    return(__mm_exactbin_w(x, w, g, dir, include))
}

real colvector __mm_exactbin(
    real colvector x,    // data (assumed sorted)
    real colvector g,    // grid (assumed sorted)
    real scalar dir,     // 0 left inclusive (default),  else right inclusive
    real scalar include) // include data outside grid in first/last bin
{
    real scalar     i, i1, j, gj
    real colvector  c

    j = rows(g)
    if (j<2) _error(3200)
    i = rows(x)
    if (i<1) return(J(j-1, 1, 0))
    if (include==0) {
        if (x[1]<g[1] | x[i]>g[j]) _error("data out of grid range")
    }
    if (j<3) return(i)
    j--
    c = J(j, 1, 0)
    if (dir==0) {
        for (; j; j--) {
            gj = g[j]
            i1 = i
            for (;i;i--) {
                if (x[i]<gj) break
            }
            if (i1<=i) continue // no data in current bin
            c[j] = i1 - i
        }
    }
    else {
        for (; j; j--) {
            gj = g[j]
            i1 = i
            for (;i;i--) {
                if (x[i]<=gj) break
            }
            if (i1<=i) continue // no data in current bin
            c[j] = i1 - i
        }
    }
    if (i) c[1] = c[1] + i // data below grid
    return(c)
}

real colvector __mm_exactbin_w(
    real colvector x,    // data (assumed sorted)
    real colvector w,    // weights
    real colvector g,    // grid (assumed sorted)
    real scalar dir,     // 0 left inclusive (default),  else right inclusive
    real scalar include) // include data outside grid in first/last bin
{
    real scalar     i, i1, j, gj
    real colvector  c
    
    j = rows(g)
    if (j<2) _error(3200)
    i = rows(x)
    if (i<1) return(J(j-1, 1, 0))
    if (include==0) {
        if (x[1]<g[1] | x[i]>g[j]) _error("data out of grid range")
    }
    if (j<3) return(sum(w))
    j--
    c = J(j, 1, 0)
    if (dir==0) {
        for (; j; j--) {
            gj = g[j]
            i1 = i
            for (;i;i--) {
                if (x[i]<gj) break
            }
            if (i1<=i) continue // no data in current bin
            c[j] = sum(w[|i+1 \ i1|])
        }
    }
    else {
        for (; j; j--) {
            gj = g[j]
            i1 = i
            for (;i;i--) {
                if (x[i]<=gj) break
            }
            if (i1<=i) continue // no data in current bin
            c[j] = sum(w[|i+1 \ i1|])
        }
    }
    if (i) c[1] = c[1] + sum(w[|1 \ i|]) // data below grid
    return(c)
}

real colvector mm_fastexactbin(
    real colvector x,      // data points
    real colvector w,      // weights
    real colvector g,      // grid (assumed sorted)
    | real scalar dir,     // 0 left inclusive (default),  else right inclusive
      real scalar include) // include data outside grid in first/last bin
{
    real scalar     g1, i, j, n, d, xi
    real colvector  c
    real rowvector  minmax
    pointer scalar  wi

    if (args()<4) dir = 0
    if (args()<5) include = 0
    n = rows(g) - 1
    i = rows(x)
    if (i<1) return(J(n, 1, 0))
    if (include==0) {
        minmax = minmax(x)
        if (minmax[1]<g[1] | minmax[2]>g[n+1]) _error("data out of grid range")
    }
    if (n<2) return(mm_nobs(x,w))
    wi = (rows(w)==1 ? &1 : (rows(w)!=i ? _error(3200) : &i))
    g1 = g[1]
    d = (g[n+1] - g1)/n
    c = J(n, 1, 0)
    if (dir==0) {
        for (; i; i--) {
            xi = x[i]
            j = floor((xi - g1) / d) + 1
            // data below grid
            if (j<1) {
                c[1] = c[1] +  w[*wi]
                continue
            }
            // data above grid
            if (j>n) {
                c[n] = c[n] +  w[*wi]
                continue
            }
            // roundoff error 1: j erroneously rounded up
            if (xi<g[j]) {
                // not sure whether this can happen at all...
                if (j>1) c[j-1] = c[j-1] + w[*wi]
                else     c[j]   = c[j] + w[*wi]
                continue
            }
            // roundoff error 2: j erroneously rounded down
            if (xi>=g[j+1]) {
                if (j<n) c[j+1] = c[j+1] + w[*wi]
                else     c[j]   = c[j] + w[*wi]
                continue
            }
            // no roundoff error
            c[j] = c[j] + w[*wi]
        }
    }
    else {
        for (; i; i--) {
            xi = x[i]
            j = ceil((xi - g1) / d)
            // data below grid
            if (j<1) {
                c[1] = c[1] +  w[*wi]
                continue
            }
            // data above grid
            if (j>n) {
                c[n] = c[n] +  w[*wi]
                continue
            }
            // roundoff error 1: j erroneously rounded up
            if (xi<=g[j]) {
                if (j>1) c[j-1] = c[j-1] + w[*wi]
                else     c[j]   = c[j] + w[*wi]
                continue
            }
            // roundoff error 2: j erroneously rounded down
            if (xi>g[j+1]) {
                // not sure whether this can happen at all...
                if (j<n) c[j+1] = c[j+1] + w[*wi]
                else     c[j]   = c[j] + w[*wi]
                continue
            }
            // no roundoff error
            c[j] = c[j] + w[*wi]
        }
    }
    return(c)
}

end
