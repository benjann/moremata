*! version 1.0.0  25jul2021  Ben Jann
version 11.2

mata:

real colvector mm_wbal(real matrix X, real colvector w, 
    real matrix X0, real colvector w0, | real scalar nowarn)
{
    class mm_ebalance scalar S
    
    if (args()<5) nowarn = 0
    S.nowarn(nowarn)
    S.trace("none")
    S.data(X, w, X0, w0)
    return(S.wbal())
}

end

