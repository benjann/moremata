*! version 1.0.7  03aug2007  Ben Jann
version 9.2
mata:

real matrix mm_ecdf(real matrix X, | real colvector w, real scalar mid)
{
    real scalar W

    if (args()<2) w = 1
    if (args()<3) mid = 0

    return(mm_ranks(X, w, 3, mid, 1))
}

end
