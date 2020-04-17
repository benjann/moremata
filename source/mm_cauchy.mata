version 9.2
// CFBaum 3aug2007
// minor changes by Ben Jann, 09aug2007
// from Wikipedia
// http://en.wikipedia.org/wiki/Cauchy_distribution
mata:
// density function of the Cauchy-Lorentz distribution
real matrix mm_cauchyden(real matrix x, real scalar x0, real scalar gamma)
{
    real matrix cauchyden

    if (gamma <= 0) {
        return(x*.)
    }
    cauchyden = 1 :/ ( pi() * gamma :* ( 1 :+ ( (x :- x0) :/ gamma ) :^2) )
    return(cauchyden)
}

// cumulative distribution function of the Cauchy-Lorentz distribution
real matrix mm_cauchy(real matrix x, real scalar x0, real scalar gamma)
{
    real matrix cauchy

    if (gamma <= 0) {
        return(x*.)
    }
    cauchy = 1 / pi() :* atan( (x :- x0) :/ gamma ) :+ 0.5
    return(cauchy)
}

// tail probability of the Cauchy-Lorentz distribution
// note: cauchy(x,0,1) is the standard Cauchy distribution
//       with cauchytail(x,0,1) = ttail(1,x)
real matrix mm_cauchytail(real matrix x, real scalar x0, real scalar gamma)
{
    return(1 :- mm_cauchy(x, x0, gamma))
}

// inverse cumulative distribution function of the Cauchy-Lorentz distribution
real matrix mm_invcauchy(real matrix p, real scalar x0, real scalar gamma)
{
    real matrix invcauchy

    if ( gamma <= 0 | min(p) < 0 | max(p) > 1 ) {
        return(p*.)
    }
    invcauchy = x0 :+ gamma :* tan( pi() :* ( p :- 0.5 ) ) :+ ln(p:>=0:&p:<=1)
    return(invcauchy)
}
end
