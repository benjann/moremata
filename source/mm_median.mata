*! version 1.1.3  Ben Jann  14jan2022
version 9.2
mata:

real rowvector mm_median(real matrix X, | real colvector w, real scalar def, 
    real scalar fw, real scalar wd)
{
    if (args()<2) w = 1
    if (args()<3) def = 2
    if (args()<4) fw = 0
    return(mm_quantile(X, w, .5, def, fw, wd))
}

real scalar _mm_median(real colvector X, | real colvector w, real scalar def, 
    real scalar fw, real scalar wd)
{   // X assumed sorted and non-missing
    // w assumed non-missing and non-negative
    if (args()<2) w = 1
    if (args()<3) def = 2
    if (args()<4) fw = 0
    return(_mm_quantile(X, w, .5, def, fw, wd))
}

end
