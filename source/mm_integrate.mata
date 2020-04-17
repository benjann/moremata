*! version 1.0.0  Ben Jann  26jan2014

version 9.2

local opts "o1, o2, o3, o4, o5, o6, o7, o8, o9, o10"

mata:

// Simpson's rule (quadratic interpolation)
real colvector mm_integrate_sr(
    pointer(real scalar function) scalar f, // function to be integrated
    real scalar a,      // lower limit
    real scalar b,      // upper limit
    real scalar n,      // integration intervals
    | real scalar el,   // element wise evaluation
    `opts')             // additional args to pass on to f
{
    transmorphic    fs            // setup for f
    real scalar     d, i
    real colvector  x, w
    
    // check options
    if (missing(n))       _error(3300, "n: missing not allowed")
    if ((n<2) | mod(n,2)) _error(3300, "n: must be a positive multiple of two")
    if (missing(a))       _error(3300, "a: missing not allowed")
    if (missing(b))       _error(3300, "b: missing not allowed")
    if (a>=b)             _error(3300, "a must be strictly smaller than b")
    
    // prepare function call
    fs = mm_callf_setup(f, args()-5, `opts')
    
    // integration setup
    d = (b-a)/n
    x = rangen(a, b, n+1)
    w = 1 \ colshape((J(n/2,1,4), J(n/2,1,2)), 1)
    w[n+1] = 1
    
    // pass column vector
    if (el==0 | args()<5) {
        return(d/3 * quadcolsum(mm_callf(fs, x) :* w))
    }
    // element wise evaluation
    for (i=1; i<=(n+1); i++) {
        x[i] = mm_callf(fs, x[i])
    }
    return(d/3 * quadcolsum(x :* w))
}

// Simpson's 3/8 rule (cubic interpolation)
real colvector mm_integrate_38(
    pointer(real scalar function) scalar f, // function to be integrated
    real scalar a,      // lower limit
    real scalar b,      // upper limit
    real scalar n,      // integration intervals
    | real scalar el,   // element wise evaluation
    `opts')             // additional args to pass on to f
{
    transmorphic    fs            // setup for f
    real scalar     d, i
    real colvector  x, w
    
    // check options
    if (missing(n))       _error(3300, "n: missing not allowed")
    if ((n<3) | mod(n,3)) _error(3300, "n: must be a positive multiple of three")
    if (missing(a))       _error(3300, "a: missing not allowed")
    if (missing(b))       _error(3300, "b: missing not allowed")
    if (a>=b)             _error(3300, "a must be strictly smaller than b")
    
    // prepare function call
    fs = mm_callf_setup(f, args()-5, `opts')
    
    // integration setup
    d = (b-a)/n
    x = rangen(a, b, n+1)
    w = 1 \ colshape((J(n/3,1,3), J(n/3,1,3), J(n/3,1,2)), 1)
    w[n+1] = 1
    
    // pass column vector
    if (el==0 | args()<5) {
        return(3*d/8 * quadcolsum(mm_callf(fs, x) :* w))
    }
    // element wise evaluation
    for (i=1; i<=(n+1); i++) {
        x[i] = mm_callf(fs, x[i])
    }
    return(3*d/8 * quadcolsum(x :* w))
}

end
