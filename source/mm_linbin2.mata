*! version 1.0.0  27may2025  Ben Jann
version 9.2
mata:

struct mm_linbin2_struct {
    real matrix    y  // aggregated data
    real colvector x  // grid points
    real colvector w  // aggregated weights
    real scalar    n  // grid size
    real scalar    d  // step size
}

struct mm_linbin2_struct scalar mm_linbin2(
    real matrix    y, // data to be aggregated; can be empty
    real colvector x, // running variable
    real colvector w, // weights; can be 1
    real scalar    n, // grid size
  | real vector    e) // outer padding of grid or (min,max); default 0
{
    real scalar flip
    struct mm_linbin2_struct scalar S
    
    if (args()<5) e = 0
    flip = 0
    // initialize
    S.n = trunc(n)
    if (length(e)==2) S.x = rangen(e[1], e[2], S.n) // e = (min,max)
    else              S.x = _mm_linbin2_grid(x, S.n, e)
    S.d = (S.x[S.n] - S.x[1]) / (S.n - 1)
    if (S.n<1 |!rows(x)) {      // empty grid or empty data
        S.y = J(S.n, cols(y), 0)
        S.w = J(S.n, 1, 0)
        return(S)
    }
    if (S.n<2) {                // one grid point only
        S.y = mean(y, w)
        S.w = mm_nobs(x, w)
        return(S)
    }
    if (S.d==0 | S.d>=.) {      // invalid grid
        S.y     = J(S.n, cols(y), 0)
        S.w     = J(S.n, 1, 0)
        S.y[1,] = mean(y, w)    // assign all data to first point
        S.w[1]  = mm_nobs(x, w)
        return(S)
    }
    if (S.d<0) {                // reverse grid
        flip = 1
        S.x = S.x[S.n::1]
        S.d = (S.x[S.n] - S.x[1]) / (S.n - 1)
    }
    // aggregate data
    if (rows(w)==1) {
        _mm_linbin2_1(S, y, x)
        if (w!=1) S.w = S.w * w
    }
    else _mm_linbin2_1w(S, y, x, w)
    _mm_linbin2_avg(S)
    if (flip) {
        S.y = S.y[S.n::1,]
        S.x = S.x[S.n::1]
        S.w = S.w[S.n::1]
    }
    return(S)
}

real colvector _mm_linbin2_grid(real colvector x, real scalar n, real scalar e)
{
    real rowvector minmax
    
    minmax = minmax(x)
    return(rangen(minmax[1] - e, minmax[2] + e, n))
}

void _mm_linbin2_avg(struct mm_linbin2_struct scalar S)
{   // compute average within each bin (unless weight is zero)
    real scalar  i
    
    for (i=S.n; i; i--) {
        if (S.w[i]) S.y[i,] = S.y[i,] :/ S.w[i]
    }
}

void _mm_linbin2_1(struct mm_linbin2_struct scalar S, real matrix y,
    real colvector x)
{
    real scalar    n, i, c, z, j
    real colvector Z, J, Sy, Sw

    n  = S.n
    c  = cols(y)
    Sy = J(n, c, 0)
    Sw = J(n, 1, 0)
    Z  = (x :- S.x[1]) / S.d
    J  = floor(Z) :+ 1
    Z  = J - Z
    for (i=rows(x); i; i--) {
        j = J[i]
        if      (j<1)  j = 1
        else if (j>=n) j = n
        else {
            z = Z[i]
            Sw[j] = Sw[j] + z
            if (c) Sy[j,] = Sy[j,] + z * y[i,]
            j++
            z = 1 - z
            Sw[j] = Sw[j] + z
            if (c) Sy[j,] = Sy[j,] + z * y[i,]
            continue
        }
        Sw[j] = Sw[j] + 1
        if (c) Sy[j,] = Sy[j,] + y[i,]
    }
    S.y = Sy; S.w = Sw 
}

void _mm_linbin2_1w(struct mm_linbin2_struct scalar S, real matrix y,
    real colvector x, real colvector w)
{
    real scalar    n, i, c, z, j, f, wi
    real colvector Z, J, Sy, Sw
    
    n  = S.n
    c  = cols(y)
    Sy = J(n, c, 0)
    Sw = J(n, 1, 0)
    Z  = (x :- S.x[1]) / S.d
    J  = floor(Z) :+ 1
    Z  = J - Z
    for (i=rows(x); i; i--) {
        j = J[i]; wi = w[i]
        if      (j<1)  j = 1
        else if (j>=n) j = n
        else {
            z = Z[i]
            f = wi * z
            Sw[j] = Sw[j] + f
            if (c) Sy[j,] = Sy[j,] + f * y[i,]
            j++
            f = wi * (1 - z)
            Sw[j] = Sw[j] + f
            if (c) Sy[j,] = Sy[j,] + f * y[i,]
            continue
        }
        Sw[j] = Sw[j] + wi
        if (c) Sy[j,] = Sy[j,] + wi * y[i,]
    }
    S.y = Sy; S.w = Sw 
}

end
