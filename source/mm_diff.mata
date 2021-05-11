*! version 1.0.1  Ben Jann  11may2021
version 9.2
mata:

numeric vector mm_diff(numeric vector x, | real scalar lag0)
{
    real scalar n, lag
    
    if (lag0>=.) lag =  1
    else         lag = abs(trunc(lag0))
    n = length(x)
    if (n<=lag) {
        if (cols(x)!=1) return(J(1,0,missingof(x)))
        return(J(0,1,missingof(x)))
    }
    return(x[|1+lag \ .|] - x[|1 \ n-lag|])
}

numeric matrix mm_rowdiff(numeric matrix x, | real scalar lag0)
{
    real scalar n, lag
    
    if (lag0>=.) lag =  1
    else         lag = abs(trunc(lag0))
    n = cols(x)
    if (n<=lag)     return(J(rows(x),0,missingof(x)))
    if (rows(x)==0) return(J(0, n-lag, missingof(x)))
    return(x[|1,1+lag \ .,.|] - x[|1,1 \ .,n-lag|])
}

numeric matrix mm_coldiff(numeric matrix x, | real scalar lag0)
{
    real scalar n, lag
    
    if (lag0>=.) lag =  1
    else         lag = abs(trunc(lag0))
    n = rows(x)
    if (n<=lag)     return(J(0,cols(x),missingof(x)))
    if (cols(x)==0) return(J(n-lag, 0, missingof(x)))
    return(x[|1+lag,1 \ .,.|] - x[|1,1 \ n-lag,.|])
}

end
