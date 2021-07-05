*! version 1.0.0  Ben Jann  04jul2021
version 9.2
mata:

real matrix mm_hdq(real matrix X, | real colvector w, real matrix P,
    real scalar fw)
{
    if (args()<2) w = 1
    if (args()<3) P = (0, .25, .50, .75, 1)'
    if (args()<4) fw = 0
    return(mm_quantile(X, w, P, 10, fw))
}

real matrix _mm_hdq(real colvector X, | real colvector w, real matrix P,
    real scalar fw)
{   // X assumed sorted and non-missing
    // w assumed non-negative and non-missing
    if (args()<2) w = 1
    if (args()<3) P = (0, .25, .50, .75, 1)'
    if (args()<4) fw = 0
    return(_mm_quantile(X, w, P, 10, fw))
}

real rowvector mm_hdmed(real matrix X, | real colvector w, real scalar fw)
{
    if (args()<2) w = 1
    if (args()<3) fw = 0
    return(mm_quantile(X, w, .5, 10, fw))
}

real scalar _mm_hdmed(real colvector X, | real colvector w, real scalar fw)
{   // X assumed sorted and non-missing
    // w assumed non-negative and non-missing
    if (args()<2) w = 1
    if (args()<3) fw = 0
    return(_mm_quantile(X, w, .5, 10, fw))
}

real matrix mm_hdiqr(real matrix X, | real colvector w, real scalar fw)
{
    real matrix q

    if (args()<2) w = 1
    if (args()<3) fw = 0
    q = mm_quantile(X, w, (.25 \ .75), 10, fw)
    return(q[2,] - q[1,])
}

real scalar _mm_hdiqr(real colvector X, | real colvector w, real scalar fw)
{   // X assumed sorted and non-missing
    // w assumed non-negative and non-missing
    real matrix q

    if (args()<2) w = 1
    if (args()<3) fw = 0
    q = _mm_quantile(X, w, (.25 \ .75), 10, fw)
    return(q[2] - q[1])
}

end
