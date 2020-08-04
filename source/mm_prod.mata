*! version 1.0.0  Ben Jann  04aug2020
*! based on a Statalist post by Daniel Klein, see:
*! www.statalist.org/forums/forum/general-stata-discussion/general/1330558-product-of-row-elements

version 10.1
mata:

real scalar mm_prod(real matrix x)
{
    if (length(x)==0) return(.)
    return(exp(sum(ln(abs(x)), 1)) * (1:-2*mod(sum(x:<0), 2)))
}

real colvector mm_rowprod(real matrix x)
{
    if (length(x)==0) return(J(rows(x), 1, .))
    return(exp(rowsum(ln(abs(x)), 1)) :* (1:-2*mod(rowsum(x:<0), 2)))
}

real rowvector mm_colprod(real matrix x) return(mm_rowprod(x')')

end
