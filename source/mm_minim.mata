*! version 1.0.0  Ben Jann  05aug2020
*! translation of function Brent_fmin() from file optimize.c included in the
*! source of R 4.0.2

/*
Quoted from optimize.c:
    an approximation  x  to the point where  f  attains a minimum  on
the interval  (ax,bx)  is determined.

INPUT..

ax    left endpoint of initial interval
bx    right endpoint of initial interval
f     function which evaluates  f(x, info)  for any  x
      in the interval  (ax,bx)
tol   desired length of the interval of uncertainty of the final
      result ( >= 0.)

OUTPUT..

fmin  abcissa approximating the point where  f  attains a minimum

    The method used is a combination of  golden  section  search  and
successive parabolic interpolation.  convergence is never much slower
than  that  for  a  Fibonacci search.  If  f  has a continuous second
derivative which is positive at the minimum (which is not  at  ax  or
bx),  then  convergence  is  superlinear, and usually of the order of
about  1.324....
    The function  f  is never evaluated at two points closer together
than  eps*abs(fmin)+(tol/3), where eps is  approximately  the  square
root  of  the  relative  machine  precision.   if   f   is a unimodal
function and the computed values of   f   are  always  unimodal  when
separated  by  at least  eps*abs(x)+(tol/3), then  fmin  approximates
the abcissa of the global minimum of  f  on the interval  ax,bx  with
an error less than  3*eps*abs(fmin)+tol.  if   f   is  not  unimodal,
then fmin may approximate a local, but perhaps non-global, minimum to
the same accuracy.
    This function subprogram is a slightly modified  version  of  the
Algol  60 procedure  localmin  given in Richard Brent, Algorithms for
Minimization without Derivatives, Prentice-Hall, Inc. (1973).
*/

version 9.2

local opts "o1, o2, o3, o4, o5, o6, o7, o8, o9, o10"

mata:

real scalar mm_minim(
    pointer(real scalar function) scalar f, // address of function to be minimized
    real scalar ax,   // minimum will be sought for within [ax,bx]
    real scalar bx,
  | real scalar tol0, // acceptable tolerance (length of uncertainty interval)
    `opts')           // additional args to pass on to f
{
    transmorphic  fs            // setup for f
    real scalar   a, b, c, d, e, p, q, r, u, v, w, x
    real scalar   t2, fu, fv, fw, fx, xm, eps, tol1, tol3, tol
    
    if (tol0<=0) _error(3300)
    tol = (tol0>=. ? epsilon(1)^0.25 : tol0) // default tolerance
    fs = mm_callf_setup(f, args()-4, `opts') // prepare function call
    c = (3 - sqrt(5)) * .5                   // squared inverse of the golden ratio
    eps = epsilon(1)
    tol1 = eps + 1 // the smallest 1.000... > 1
    eps = sqrt(eps)
    a = ax; b = bx
    v = w = x = a + c * (b - a)
    d = e = 0
    fv = fw = fx = mm_callf(fs, x)
    tol3 = tol / 3
    while (1) {
        xm = (a + b) * .5
        tol1 = eps * abs(x) + tol3
        t2 = tol1 * 2
        /* check stopping criterion */
        if (abs(x - xm) <= t2 - (b - a) * .5) break
        p = q = r = 0
        /* fit parabola */
        if (abs(e) > tol1) {
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2
            if (q > 0) p = -p
            else       q = -q
            r = e
            e = d
        }
        /* a golden-section step */
        if (abs(p) >= abs(q * .5 * r) | p <= q * (a - x) | p >= q * (b - x)) {
            if (x < xm) e = b - x
            else        e = a - x
            d = c * e
        }
        /* a parabolic-interpolation step */
        else {
            d = p / q
            u = x + d
            /* f must not be evaluated too close to ax or bx */
            if (u - a < t2 | b - u < t2) {
                d = tol1
                if (x >= xm) d = -d
            }
        }
        /* f must not be evaluated too close to x */
        if (abs(d) >= tol1) u = x + d
        else if (d > 0)     u = x + tol1
        else                u = x - tol1
        fu = mm_callf(fs, u)
        /*  update  a, b, v, w, and x */
        if (fu <= fx) {
            if (u < x) b = x
            else       a = x
            v = w;    w = x;   x = u
            fv = fw; fw = fx; fx = fu
        }
        else {
            if (u < x) a = u
            else       b = u
            if (fu <= fw | w == x) {
                v = w; fv = fw
                w = u; fw = fu
            } 
            else if (fu <= fv | v == x | v == w) {
                v = u; fv = fu
            }
        }
    }
    return(x)
}

end
