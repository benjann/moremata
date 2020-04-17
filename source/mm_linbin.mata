*! version 1.0.8  03aug2007  Ben Jann
version 9.0
mata:

real colvector mm_linbin(
 real colvector x,  // data points
 real colvector w,  // weights
 real colvector g)  // grid points
{
    real colvector  c, p, ng, nx
    real scalar     i, j, d, ww


    ng = rows(g)
    nx = rows(x)
    ww = rows(w)!=1
    if ((ww & rows(w)!=nx) | ng<1) _error(3200)
    if (nx<1) return(J(ng, 1, 0))
    if (ng<2) return(ww ? colsum(w) : w*nx)
    p = order((x,(1::nx)),(1,2)) //stable sort order
    c = J(ng, 1, 0)
    for (i=1; i<=nx; i++) { // data below grid
        if (x[p[i]]>=g[1])  break
        c[1] = c[1] + (ww ? w[p[i]] : w)
    }
    j = 2
    for (; i<=nx; i++) {
        if (g[ng]<x[p[i]]) break
        while (g[j]<x[p[i]]) j++
        d = (ww ? w[p[i]] : w) / (g[j] - g[j-1])
        c[j-1] = c[j-1] + (g[j] - x[p[i]]) * d
        c[j]   = c[j] + (x[p[i]] - g[j-1]) * d
    }
    if (i<=nx) { // data above grid
        c[ng] = c[ng] + (ww ? colsum(w[p[|i \ .|]]) : w*(nx-i+1))
    }
    return(c)
}

end
