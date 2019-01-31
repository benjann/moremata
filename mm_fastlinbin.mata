*! version 1.0.2  03aug2007  Ben Jann
version 9.1
mata:

real colvector mm_fastlinbin(
 real colvector x,  // data points
 real colvector w,  // weights
 real colvector g)  // grid points
{
    real colvector  c
    real scalar     g1, i, j, d, z, ng, nx
    pointer scalar  wi

    ng = rows(g)
    nx = rows(x)
    if (ng<1) _error(3200)
    if (nx<1) return(J(ng, 1, 0))
    if (ng<2) return(mm_nobs(x,w))
    wi = (rows(w)==1 ? &1 : (rows(w)!=nx ? _error(3200) : &i))
    g1 = g[1]
    d = (g[ng] - g1)/(ng-1)
    c = J(ng, 1, 0)
    for (i=1; i<=rows(x); i++) {
        z = (x[i] - g1) / d
        if (z<=0) {
            c[1] = c[1] +  w[*wi]
            continue
        }
        if (z>=(ng-1)) {
            c[ng] = c[ng] +  w[*wi]
            continue
        }
        j = trunc(z) + 1
        c[j] = c[j] + w[*wi] * (j-z)
        j++
        c[j] = c[j] + w[*wi] * (2-j+z)
    }
    return(c)
}
end
