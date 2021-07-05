*! version 1.1.0  Ben Jann  13jul2020
version 9.2
mata:

real matrix mm_iqrange(real matrix X, | real colvector w, real scalar def, 
    real scalar fw)
{
    real matrix q

    if (args()<2) w = 1
    if (args()<3) def = 2
    if (args()<4) fw = 0
    q = mm_quantile(X, w, (.25 \ .75), def, fw)
    return(q[2,] - q[1,])
}

real scalar _mm_iqrange(real colvector X, | real colvector w, real scalar def, 
    real scalar fw)
{   // X assumed sorted and non-missing
    // w assumed non-missing and non-negative
    real matrix q

    if (args()<2) w = 1
    if (args()<3) def = 2
    if (args()<4) fw = 0
    q = _mm_quantile(X, w, (.25 \ .75), def, fw)
    return(q[2] - q[1])
}

end
