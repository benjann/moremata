*! version 1.0.1  16jun2025  Ben Jann
version 14.2
mata:

real vector mm_loclin(real colvector Y, real colvector X, real colvector w,
    real vector at, real vector bw, | string scalar kernel)
{
    real scalar    i, varbw
    string scalar  kern
    real vector    fit
    pointer scalar k
    
    kern = _mm_unabkern(kernel) // default is "epan2"
    varbw = length(bw)!=1
    fit = J(rows(at), cols(at), .)
    if (kern=="gaussian") {
        for (i=length(fit); i; i--) {
            fit[i] = _mm_loclin_gaussian(Y, X, w, at[i], varbw ? bw[i] : bw)
        }
        return(fit)
    }
    k = _mm_findkern(kern)
    for (i=length(fit); i; i--) {
        fit[i] = _mm_loclin(Y, X, w, at[i], varbw ? bw[i] : bw, k)
    }
    return(fit)
}

real scalar _mm_loclin_gaussian(real colvector Y, real colvector X,
    real colvector w, real scalar at, real scalar bw)
{
    real colvector Z
    
    Z = X :- at
    return(mm_lsfit(Y, Z, w :* mm_kern_gaussian(Z/bw))[2])
}

real scalar _mm_loclin(real colvector Y, real colvector X, real colvector w,
    real scalar at, real scalar bw, pointer scalar kern)
{
    real scalar    n
    real colvector W, Z, p
    
    Z = X :- at
    W = (*kern)(Z/bw)
    p = selectindex(W)
    n = length(p)
    if (!n)         return(.) // no data within window
    if (n==rows(Y)) return(mm_lsfit(Y, Z, W:*w)[2])
                    return(mm_lsfit(Y[p], Z[p], (W:*w)[p])[2])
}

real scalar mm_loclin_bw(real colvector Y, real colvector X,
    | real colvector w, string scalar kernel, real scalar fw)
{
    real scalar    h, c
    string scalar  kern
    pointer scalar k
    
    if (args()<3) w  = 1
    if (args()<5) fw = 0
    // compute canonical bandwidth of kernel; using mm_kdel0(kernel) would be
    // more precise and more efficient, but the goal is to reproduce lpoly's
    // bandwidth as closely as possible; thanks to Jeff Pitblado from Stata
    // Corp for pointing me to the existence of (undocumented) function
    // _lpoly_CnupK()
    kern = _mm_unabkern(kernel) // default is "epan2"
    if (kern=="epanechnikov")   {; h = sqrt(5); k = &_lpoly_epanechnikov(); }
    else if (kern=="epan2")     {; h = 1.25;    k = &_lpoly_epan2(); }
    else if (kern=="biweight")  {; h = 1.25;    k = &_lpoly_biweight(); }
    else if (kern=="cosine")    {; h = 1.25;    k = &_lpoly_cosine(); }
    else if (kern=="gaussian")  {; h = 5;       k = &_lpoly_gaussian(); }
    else if (kern=="parzen")    {; h = 1.25;    k = &_lpoly_parzen(); }
    else if (kern=="rectangle") {; h = 1.25;    k = &_lpoly_rectangle(); }
    else if (kern=="triangle")  {; h = 1.25;    k = &_lpoly_triangle(); }
    else _error(3300) // (cannot be reached)
    c = _lpoly_CnupK(-h, h, 2*h/100, 1, 0, k)
    // compute ROT bandwiath
    return(_mm_loclin_bw(Y, X, w, fw, 0) * c)
}

real scalar mm_loclin_bw2(real colvector Y, real colvector X,
    | real colvector w, string scalar kernel, real scalar fw)
{
    if (args()<3) w  = 1
    if (args()<5) fw = 0
    return(_mm_loclin_bw(Y, X, w, fw, 1) * mm_kdel0(kernel))
}

real scalar _mm_loclin_bw(real colvector Y, real colvector X, real colvector w,
    real scalar fw, real scalar alt)
{
    real scalar    ll, ul, s2, M2, n
    real rowvector minmax
    real colvector b
    transmorphic   S
    
    minmax = minmax(X)
    ll = minmax[1] + .05*(minmax[2] - minmax[1])
    ul = minmax[2] - .05*(minmax[2] - minmax[1])
    S  = mm_ls(Y, (X,X:^2,X:^3,X:^4), w)
    b  = mm_ls_b(S)
    s2 = mm_ls_rss(S)/mm_ls_N(S)
    n  = fw ? (rows(w)==1 ? w*rows(X) : sum(w)) : rows(X)
    // Stata's lpoly takes the integral of the 2nd derivative of the global
    // polynomial fit without considering how X is distributed (i.e. assuming
    // an even distribution); however, according to the definition of the ROT
    // bandwidth by Fan and Gijbels (1996, page 111) the integral should be
    // taken across the distribution of X; this is what is done if alt!=0; else
    // the same approach as in Stata's lpoly is used
    if (alt) M2 = mean(_mm_loclin_bw_d2(X, b) :* (X:>=ll :& X:<=ul), w)
    else     M2 = mm_integrate_sr(&_mm_loclin_bw_d2(), ll, ul, 1000, 0, b)
    return((s2 * (ul-ll) / (n * M2))^0.2)
}

real colvector _mm_loclin_bw_d2(real colvector X, real vector b)
    return((2*b[2] :+ 6*b[3]*X :+ 12*b[4]*X:^2):^2)

end
