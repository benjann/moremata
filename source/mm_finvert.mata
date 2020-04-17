*! version 1.0.1, Ben Jann, 16may2014

version 9.2

mata:

real vector mm_finvert(
    y,
    f,     //     nr:      brent:
    o1,    //     fd       lo
    | o2,  //     x0       up
      tol,
      maxit,
      opt)
{
    if (args()<5) tol = 0
    if (args()<6) maxit = 1000
    if (ispointer(o1)) {
        if (args()<4) o2 = y   // set initial guess to y
        if (args()<7)
            return(_mm_finvert_nr(y, f, o1, o2, tol, maxit))
        else
            return(_mm_finvert_nr(y, f, o1, o2, tol, maxit, opt))
    }
    if (args()==3)
        _error(3001, "expected 4 to 6 arguments but received 3")
    if (args()<7)
        return(_mm_finvert_brent(y, f, o1, o2, tol, maxit))
    else 
        return(_mm_finvert_brent(y, f, o1, o2, tol, maxit, opt))
}

real vector _mm_finvert_nr(
    real vector y,
    pointer(real scalar function) scalar f,
    pointer(real scalar function) scalar df,
    real vector x,
    real scalar tol,
    real scalar maxit,
    | opt)
{
    real scalar     i, I, resi, rc
    real vector     res
    pointer scalar  ix

    I = length(y)
    ix = (length(x)==1 ? &1 : (length(x)<I ? _error(3200) : &i))

    res = J(rows(y), cols(y), .)
    for (i=1; i<=I; i++) {
        if (args()<7)
            rc = mm_nrroot(resi=x[*ix], &_mm_finvert_nr_fdf(),
                tol, maxit, y[i], f, df)
        else
            rc = mm_nrroot(resi=x[*ix], &_mm_finvert_nr_fdf(),
                tol, maxit, y[i], f, df, opt)
        if (rc)
         _error(3360, "failure to converge for element " + strofreal(i))
        res[i] = resi
    }
    return(res)
}

real rowvector _mm_finvert_nr_fdf(
    real scalar x,
    real scalar y,
    pointer(real scalar function) scalar f,
    pointer(real scalar function) scalar df,
    | opt)
{
    if (args()<5)
        return( (*f)(x) - y, (*df)(x) )
    else
        return( (*f)(x, opt) - y, (*df)(x, opt) )
}

real vector _mm_finvert_brent(
    real vector y,
    pointer(real scalar function) scalar f,
    real vector lo,
    real vector up,
    real scalar tol,
    real scalar maxit,
    | opt)
{
    real scalar     i, I, resi, rc
    real vector     res
    pointer scalar  il, iu

    I = length(y)
    il = (length(lo)==1 ? &1 : (length(lo)<I ? _error(3200) : &i))
    iu = (length(up)==1 ? &1 : (length(up)<I ? _error(3200) : &i))

    res = J(rows(y), cols(y), .)
    for (i=1; i<=I; i++) {
        if (args()<7)
            rc = mm_root(resi=., &_mm_finvert_brent_f(),
                lo[*il], up[*iu], tol, maxit, y[i], f)
        else
            rc = mm_root(resi=., &_mm_finvert_brent_f(),
                lo[*il], up[*iu], tol, maxit, y[i], f, opt)
        if (rc==1)
         _error(3360, "failure to converge for element " + strofreal(i))
        if (rc)
         _error(3498, "invalid search range for element " + strofreal(i))
        res[i] = resi
    }
    return(res)

}

real scalar _mm_finvert_brent_f(
    real scalar x,
    real scalar y,
    pointer(real scalar function) scalar f,
    | opt)
{
    if (args()<4)
        return( (*f)(x) - y )
    else
        return( (*f)(x, opt) - y )
}

end
