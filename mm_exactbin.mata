*! version 1.0.5  03aug2007  Ben Jann
version 9.0
mata:

real colvector mm_exactbin(
 real colvector x,   // data
 real colvector w,   // weights
 real colvector g,   // grid
 | real scalar dir,  // 0 left inclusive (default),  else right inclusive
   real scalar include) // include data outside grid in first/last bin
{
    real colvector  c, p, nx, ng
    real scalar     i, j, ww

    if (args()<4) dir = 0
    if (args()<5) include = 0
    nx = rows(x)
    ng = rows(g)
    ww = rows(w)!=1
    if ((ww & rows(w)!=nx) | rows(g)<2)  _error(3200)
    if (nx<1) return(J(ng-1, 1, 0))
    p = order((x,(1::nx)),(1,2)) //stable sort order
    if (include==0) {
        if (x[p[1]]<g[1] | x[p[nx]]>g[ng])
            _error("data out of grid range")
    }
    if (ng<3) return(ww ? colsum(w) : w*nx)
    c = J(ng-1, 1, 0)
    if (dir==0) {
        j = rows(c)
        for (i=nx; i>=1; i--) {
            while (g[j]>x[p[i]]) {
                if (j<=1) break
                j--
            }
            c[j] = c[j] + (ww ? w[p[i]] : w)
        }
    }
    else {
        j = 2
        for (i=1; i<=nx; i++) {
            while (g[j]<x[p[i]]) {
                if (j>=ng) break
                j++
            }
            c[j-1] = c[j-1] + (ww ? w[p[i]] : w)
        }
    }
    return(c)
}
end
