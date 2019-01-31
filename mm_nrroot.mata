*! version 1.0.0, Ben Jann, 07jul2006

version 9.2

local opts "o1, o2, o3, o4, o5, o6, o7, o8, o9, o10"

mata:

real scalar mm_nrroot(
 real scalar x,              // initial guess; will be replaced by solution
 pointer(function) scalar f, // the function; must return f() and f'()
 | real scalar tol,          // acceptable tolerance (default 0)
   real scalar maxit,        // maximum # of iterations (default: 1000)
   `opts')                   // additional args to pass on to f
{
    transmorphic  fs         // setup for f
    real scalar   itr
    real scalar   fx         // fx = ( f(x), f'(x) )
    real scalar   dx         // dx = fx[1] / fx[2]
    real scalar   tol_act    // actual tolerance

    if (args()<3) tol = 0      // default tolerance
    if (args()<4) maxit = 1000 // default maximum # of iterations

    fs = mm_callf_setup(f, args()-4, `opts') // prepare function call
    for (itr=1; itr<=maxit; itr++) {
        tol_act = 2e+4*epsilon(x) + tol/2
        fx = mm_callf(fs, x)
        dx = fx[1] / fx[2]
        x = x - dx
        if ( abs(dx) <= tol_act | missing(x) ) return(0)
    }
    return(1)                // convergence not reached
}

end
