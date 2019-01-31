*! version 1.0.0, Ben Jann, 07jul2006
*! translation of zeroin.c from http://www.netlib.org/c/brent.shar

* Quoted from http://www.netlib.org/c/brent.shar:
* X * Algorithm
* X *    G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
* X *    computations. M., Mir, 1980, p.180 of the Russian edition
* X *
* X *    The function makes use of the bissection procedure combined with
* X *    the linear or quadric inverse interpolation.
* X *    At every step program operates on three abscissae - a, b, and c.
* X *    b - the last and the best approximation to the root
* X *    a - the last but one approximation
* X *    c - the last but one or even earlier approximation than a that
* X *        1) |f(b)| <= |f(c)|
* X *        2) f(b) and f(c) have opposite signs, i.e. b and c confine
* X *           the root
* X *    At every step Zeroin selects one of the two new approximations, the
* X *    former being obtained by the bissection procedure and the latter
* X *    resulting in the interpolation (if a,b, and c are all different
* X *    the quadric interpolation is utilized, otherwise the linear one).
* X *    If the latter (i.e. obtained by the interpolation) point is
* X *    reasonable (i.e. lies within the current interval [b,c] not being
* X *    too close to the boundaries) it is accepted. The bissection result
* X *    is used in the other case. Therefore, the range of uncertainty is
* X *    ensured to be reduced at least by the factor 1.6

* Main changes compared to the original:
* 1. mm_root() returns an error code and stored the solution in x.
* 2. The maximum # of iteration may be set.
* 3. The tolerance argument is optional.
* 4. Additional arguments may be passed on to f.
* 5. Special action is taken if f returns missing, convergence is not
*    reached, or f(ax) and f(bx) do not have opposite signs.

version 9.2

local opts "o1, o2, o3, o4, o5, o6, o7, o8, o9, o10"

mata:

real scalar mm_root(
 transmorphic x,      // bj: will be replaced by solution
 pointer(real scalar function) scalar f,
                      // Address of the function whose zero will be sought for
 real scalar ax,      // Root will be sought for within a range [ax,bx]
 real scalar bx,      //
 | real scalar tol,   // Acceptable tolerance for the root value (default 0)
   real scalar maxit, // bj: maximum # of iterations (default: 1000)
   `opts')            // bj: additional args to pass on to f
{
    transmorphic  fs            // setup for f
    real scalar   a, b, c       // Abscissae, descr. see above
    real scalar   fa, fb, fc    // f(a), f(b), f(c)
    real scalar   prev_step     // Distance from the last but one
    real scalar   tol_act       // Actual tolerance
    real scalar   p             // Interpolation step is calcu-
    real scalar   q             // lated in the form p/q; divi-
                                // sion operations is delayed
                                // until the last moment
    real scalar   new_step      // Step at this iteration
    real scalar   t1, cb, t2
    real scalar   itr

    if (args()<5) tol = 0       // bj: set tolerance
    if (args()<6) maxit = 1000  // bj: maximum # of iterations

    fs = mm_callf_setup(f, args()-6, `opts') // bj: prepare function call

    x = .                       // bj: initialize output

    a = ax;  b = bx;  fa = mm_callf(fs, a);  fb = mm_callf(fs, b)
    c = a;  fc = fa

    if ( fa==. ) return(0)      // bj: abort if fa missing

    if ( (fa > 0 & fb > 0) |    // bj: f(ax) and f(bx) do not
         (fa < 0 & fb < 0) ) {  //     have opposite signs
        if ( abs(fa) < abs(fb) ) {
            x = a; return(2)    // bj: fa closer to zero than fb
        }
        x = b; return(3)        // bj: fb closer to zero than fa
    }

    for (itr=1; itr<=maxit; itr++) {
        if ( fb==. ) return(0)            // bj: abort if fb missing

        prev_step = b-a

        if( abs(fc) < abs(fb) ) {         // Swap data for b to be the
            a = b;  b = c;  c = a;        // best approximation
            fa = fb;  fb = fc;  fc = fa
        }

        tol_act = 2*epsilon(b) + tol/2
        new_step = (c-b)/2

        if( abs(new_step) <= tol_act | fb == 0 ) {
             x = b                        // Acceptable approx. is found
             return(0)
        }

        // Decide if the interpolation can be tried
        if( abs(prev_step) >= tol_act     // If prev_step was large enough
             & abs(fa) > abs(fb) ) {      // and was in true direction,
                                          // Interpolation may be tried
            cb = c-b
            if( a==c ) {                  // If we have only two distinct
                t1 = fb/fa                // points linear interpolation
                p = cb*t1                 // can only be applied
                q = 1.0 - t1
            }
            else {                        // Quadric inverse interpolation
                q = fa/fc;  t1 = fb/fc;  t2 = fb/fa
                p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) )
                q = (q-1.0) * (t1-1.0) * (t2-1.0)
            }
            if( p>0 )                     // p was calculated with the op-
              q = -q                      // posite sign; make p positive
            else                          // and assign possible minus to
              p = -p                      // q
            if( p < (0.75*cb*q-abs(tol_act*q)/2) // If b+p/q falls in [b,c]
               & p < abs(prev_step*q/2) ) // and isn't too large
             new_step = p/q               // it is accepted
                                          // If p/q is too large then the
                                          // bisection procedure can
                                          // reduce [b,c] range to more
                                          // extent
        }

        if( abs(new_step) < tol_act ) {   // Adjust the step to be not less
            if( new_step > 0 )            // than tolerance
             new_step = tol_act
            else
             new_step = -tol_act
        }

        a = b;  fa = fb                   // Save the previous approx.
        b = b + new_step;  fb = mm_callf(fs, b) // Do step to a new approxim.
        if( (fb > 0 & fc > 0) | (fb < 0 & fc < 0) )  {
                      // Adjust c for it to have a sign opposite to that of b
            c = a;  fc = fa
        }
    }
    x = b
    return(1)                             // bj: convergence not reached
}

end
