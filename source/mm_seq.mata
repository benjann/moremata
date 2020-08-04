*! version 1.0.0  Ben Jann  04aug2020

version 9.2
mata:

real colvector mm_seq(real scalar a, real scalar b, real scalar d)
{
    if (missing((a,b,d))) return(J(0,1,.))
    if (a<=b)             return( a :+ (0::(b-a)/abs(d)) * abs(d) )
                          return( a :- (0::(a-b)/abs(d)) * abs(d) )
}

end
