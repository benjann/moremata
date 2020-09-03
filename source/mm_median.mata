*! version 1.1.1  Ben Jann  03sep2020
version 9.2
mata:

real rowvector mm_median(real matrix X, | real colvector w)
{
    if (args()<2) w = 1
    return(mm_quantile(X, w, .5, 2))
}

real scalar _mm_median(real colvector X, | real colvector w)
{   // X assumed sorted and non-missing
    // w assumed non-missing and strictly positive
    if (args()<2) w = 1
    return(_mm_quantile(X, w, .5, 2))
}

end
