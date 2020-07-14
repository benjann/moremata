*! version 1.1.0  Ben Jann  13jul2020
version 9.2
mata:

real rowvector mm_median(real matrix X, | real colvector w, real scalar fw)
{
    if (args()<2) w = 1
    if (args()<3) fw = 0
    return(mm_quantile(X, w, .5, 2, fw))
}

real scalar _mm_median(real colvector X, | real colvector w, real scalar fw)
{   // X assumed sorted and non-missing
    // w assumed non-missing and strictly positive
    if (args()<2) w = 1
    if (args()<3) fw = 0
    return(_mm_quantile(X, w, .5, 2, fw))
}

end
