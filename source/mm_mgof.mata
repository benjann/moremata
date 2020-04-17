*! version 1.0.7  Ben Jann  12jan2008
*  mm_mgof: goodness-of-fit tests for multinomial data

version 9.2
mata:

real matrix mm_mgof(
    real colvector f,    // observed count
  | real colvector h0,   // expected distribution
    string scalar m0,    // method: "approx" (default), "mc", "rc", or "ee"
    string vector s0,    // test-stats: "x2", "lr", "cr", "mlnp", "ksmirnov"
    real scalar lambda,  // lambda for Cressie-Read statistic
    real scalar arg1,    // nfit for "approx"; reps for "mc"
    real scalar dots,    // display dots with mc/ee
    real scalar arg2)    // will be replaced; used by mgof_ee to return counter
{
    real scalar     i, n, hsum
    string scalar   m, s
    string vector   stats
    real vector     p
    pointer(real colvector) scalar       h
    pointer(real scalar function) vector s1

// parsing and defaults
    if (args()<7) dots = 0
    m = mm_strexpand(m0, ("approx", "mc", "rc", "ee"), "approx", 0)
    stats = ("x2", "lr", "cr")
    if (m!="approx") stats = (stats, "mlnp", "ksmirnov")
    if (length(s0)<1) s = "x2"
    else {
        s = J(length(s0),1,"")
        for (i=1;i<=length(s0);i++) {
            s[i] = mm_strexpand(s0[i], stats, "x2", 1)
        }
    }
    if (m!="approx" & lambda<0) _error(3498, "lambda<0 not allowed with "+m)

// check data / prepare null distribution
    if (missing(f)) _error(3351, "f has missing values")
    if (any(f:<0)) _error(3498, "f has negative values")
    if (any(f:<0)) _error(3498, "f has negative values")
    if (m=="ee" | m=="rc" | (m=="mc" & anyof(s, "mlnp"))) {
        if (f!=trunc(f)) _error(3498, "non-integer f not allowed in this context")
    }
    if (anyof(stats,"cr") & lambda<0 & anyof(f,0))
     _error(3498, "some f are 0: lambda<0 not allowed")
    if (args()<2) h0 = 1
    if (rows(f)!=rows(h0) & rows(h0)!=1) _error(3200)
    if (missing(h0)) _error(3351, "h has missing values")
    if (any(h0:<=0)) _error(3498, "some h are negative or 0")
    n = colsum(f)
    hsum = colsum(h0)
    if (n!=hsum) h = &(h0*n/hsum)
    else h = &h0

// return (0,1) if n. of obs is 0 or n. of categories < 2
    if ( n==0 | rows(f)<2 ) return(J(rows(s),1,0),J(rows(s),1,1))

// gather statistics functions
    s1 = J(rows(s),1,NULL)
    for (i=1;i<=length(s);i++) {
        if      (s[i]=="x2")        s1[i] = &_mm_mgof_x2()
        else if (s[i]=="lr")        s1[i] = &_mm_mgof_lr()
        else if (s[i]=="cr")        s1[i] = &_mm_mgof_cr()
        else if (s[i]=="mlnp")      s1[i] = &_mm_mgof_mlnp()
        else if (s[i]=="ksmirnov")  s1[i] = &_mm_mgof_ksmirnov()
        else _error(3498,s[i]+ " invalid")
    }

// approximate chi2-test
    if (m=="approx") return(_mm_mgof_approx(s1, f, *h, lambda, arg1))

// Monte Carlo test
    if (m=="mc") return(_mm_mgof_mc(s1, f, *h, lambda, arg1, dots))

// exhaustive enumeration / random composition test
    if (s=="mlnp") {
        if (m=="rc") return(_mm_mgof_rc(NULL, f, *h))
        return(_mm_mgof_ee(NULL, f, *h))
    }
    p = mm_which(s:=="mlnp")
    if (rows(p)>0) {
        p = p[1]
        p = p \ (p>1 ? (1::p-1) : J(0,1,.)) \ (p<rows(s) ? (p+1::rows(s)) : J(0,1,.))
        s1 = s1[p[| 2 \ .|]]
        p = invorder(p)
    }
    else p = (2::length(s)+1)
    if (m=="rc") return(_mm_mgof_rc(s1, f, *h, lambda, arg1)[p,])
    return(_mm_mgof_ee(s1, f, *h, lambda, anyof(s,"ksmirnov"), dots, arg2)[p,])
}

// large sample chi2-approximation test for multinomial distributions
real matrix _mm_mgof_approx(
 pointer(real scalar function) vector s, // statistics to be tested
 real colvector f,     // observed count
 real colvector h0,    // expected count
 | real scalar lambda, // lambda for Cressie-Read statistic
   real scalar nfit)   // n. of fitted parameters (default 0)
{
    real scalar                     n, i, k, nstat
    real colvector                  stat
    pointer(real colvector) scalar  h

// compute hypothetical distribution if missing
    n = colsum(f)
    k = rows(f)
    if (rows(h0)<=1) h = &(n :* J(k,1,1/k))
    else             h = &h0
// compute statistics and p-values
    nstat = length(s)
    stat = J(nstat,1,.)
    for (i=1; i<=nstat; i++) stat[i] = (*s[i])(f,*h,lambda)
    return(stat, chi2tail(rows(f) - (nfit<. ? nfit : 0) - 1,stat))
}

// Monte Carlo exact test for discrete distributions
// (based on sampling from the hypothetical distribution)
real matrix _mm_mgof_mc(
 pointer(real scalar function) vector s, // statistics to be tested
 real colvector f,     // observed count
 real colvector h0,    // expected count
 | real scalar lambda, // lambda for Cressie-Read statistic
   real scalar reps0,  // replications (default=10000)
   real scalar dots)   // display dots
{
    real scalar                     n, j, i, k, nstat, reps, uni
    real scalar                     doti, ndots
    real colvector                  fj, H, stat, p, hp
    pointer(real colvector) scalar  h

// set defaults and compute hypothetical distribution
    if (args()<5 | reps0>=.) reps = 10000
    else                     reps = reps0
    if (args()<6) dots = 0
    n = round(colsum(f)) // so that, e.g., 10 are sampled if n=9.9999...
    k = rows(f)
    if (uni = (rows(h0)==1)) {
        uni = 1
        h  = &(n :* J(k,1,1/k))
        H  = rangen(1/k,1,k)
        hp = J(k,1,1/k)
    }
    else {
        uni = mm_isconstant(h0)
        H  = mm_colrunsum(h0)
        H  = H / H[rows(H)]
        hp = (h0) / colsum(h0)
        h  = &h0
    }
// compute statistics and p-values
    nstat = length(s)
    stat = p = J(nstat,1,0)
    for (i=1; i<=nstat; i++) stat[i] = (*s[i])(f,*h, lambda, n, hp, H)
    if (reps<=0) return(stat, J(nstat,1,.))
    if (dots) {
        printf("\n{txt}Percent completed ({res}%g{txt} replications)\n",reps)
        display("{txt}0 {hline 5} 20 {hline 6} 40 {hline 6} 60 {hline 6} 80 {hline 5} 100")
        ndots = (reps>=50 ? 1 : (reps>=25 ? 2 : (reps>=10 ? 5 :
                (reps>=5 ? 10 : (reps>=2 ? 25 : 50)))))
        doti = 1
    }
    for (j=1; j<=reps; j++) {
        if (uni) fj = _mm_mgof_mc_usmpl(n, k) // = mm_srswr(n,rows(f),1)
        else     fj = _mm_mgof_mc_smpl(n, H) // = mm_upswr(n, h, 1)  (slower)
        for (i=1; i<=nstat; i++)  p[i] = p[i] +
          (stat[i] <= (*s[i])(fj,*h, lambda, n, hp, H))
        if (dots) {
            if (j*50 >= doti*ndots*reps) {
                printf(ndots*".")
                displayflush()
                doti++
            }
        }
    }
    if (dots) display("")
    return(stat, p/reps)
}
// Monte Carlo subroutine: sample from uniform multinomial
real colvector _mm_mgof_mc_usmpl(
 real scalar n, // sample size
 real scalar k) // number of categories
{
    real scalar     i
    real colvector  u, res

    u = ceil(uniform(n,1)*k)
    res = J(k,1,0)
    if (n>1.5^(k-2)) { // faster for small k relative to n
        for (i=1; i<=k; i++) res[i] = colsum(u:==i)
    }
    else {
        for (i=1;i<=n;i++) res[u[i]] = res[u[i]] + 1
    }
    return(res)
}
// Monte Carlo subroutine: sample from non-uniform multinomial
real colvector _mm_mgof_mc_smpl(
 real scalar n,    // sample size
 real colvector F) // distribution function (cumulative)
{
    real scalar     i, j
    real colvector  u, res

    if (n>2.3^(rows(F)-3)) {  // faster for small k relative to n
        u = uniform(n,1)
        res = J(rows(F),1,0)
        res[1] = sum(u:<=F[1])
        for (j=2;j<=rows(F);j++) {
            res[j] = sum(u:>F[j-1] :& u:<=F[j])
        }
        return(res)
    }
    u = sort(uniform(n,1),1)
    j=1
    res = J(rows(F),1,0)
    for (i=1;i<=n;i++) {
        while (u[i]>F[j]) j++
        res[j] = res[j]+1
    }
    return(res)
}

// Random composition exact test for discrete distributions
// -> undocumented; do not use this
//    the method is highly unreliable if the number of potential
//    compositions is large
real matrix _mm_mgof_rc(
 pointer(real scalar function) vector s0, // statistics to be tested
 real colvector f,     // observed count
 real colvector h0,    // expected count
 | real scalar lambda, // lambda for Cressie-Read statistic
   real scalar reps0)  // replications (default=10000)
{
    real scalar                           n, i, j, k, nstat, reps, lnscale
    real colvector                        fj, H, stat, statj, p, jj, se, hp
    pointer(real colvector) scalar        h
    pointer(real scalar function) vector  s

// set defaults and compute hypothetical distribution
    if (args()<5 | reps0>=.) reps = 10000
    else                     reps = reps0
    if (s0==NULL) s = &_mm_mgof_mlnp()
    else {
        if (rows(s0)!=1) s = &_mm_mgof_mlnp() \ s0
        else             s = &_mm_mgof_mlnp() , s0
    }
    n = colsum(f)
    k = rows(f)
    if (rows(h0)==1) {
        h = &(n :* J(k,1,1/k))
        H = rangen(1/k,1,k)
        hp = J(k,1,1/k)
    }
    else {
        h = &h0
        H = mm_colrunsum(*h)
        H = H / H[rows(H)]
        hp = (*h) / colsum(*h)
    }
// compute statistics and p-values
    nstat = length(s)
    stat = statj = p = jj = se = J(nstat,1,0)
    for (i=1; i<=nstat; i++) stat[i] = (*s[i])(f,*h, lambda, n, hp, H)
    if (reps<=0) return(stat, J(nstat,2,.))
    for (j=1; j<=reps; j++) {
        fj = mm_rcomposition(n, k)
        for (i=1; i<=nstat; i++) {
            statj[i] = (*s[i])(fj,*h, lambda, n, hp, H)
            if (stat[i] <= statj[i]) {
                jj[i] = jj[i] + 1
                p[i]  = p[i] + (exp(-statj[1])-p[i])/jj[i]     // mean updating
                se[i] = se[i] + (exp(-statj[1])^2-se[i])/jj[i] // mean updating
            }
        }
    }
    p = p:*jj; se = se:*jj
    // lnscale = log. of inverse sampling probability
    lnscale = lnfactorial(n+k-1) - lnfactorial(k-1) - lnfactorial(n) - ln(reps)
    se = exp( lnscale :+ ln(sqrt(reps*(1/(reps-1) * se - (p/reps):^2))) )
//  se = mm_ncompositions(n, k) * sqrt((1/(reps-1) * se - (p/reps):^2)/reps)
    p  = exp( lnscale :+ ln(p) )
//  p  = mm_ncompositions(n, k) / reps * p
    return(stat, p, se)
}

// exhaustive enumeration exact test (performs the exact
// multinomial g.o.f. test plus, optionally, exact tests
// for specified additional statistics)
// WARNING: use this test for very small samples only since
// computation time is linear to the number of compositions
// =comb(n+k-1,k-1), where n is the sample size and k is
// the number of categories.
// IMPORTANT EXCEPTION: when the null distribution is uniform,
// compution time is determined by the number of partitions
// which is magnitudes smaller than the number of compositions
real matrix _mm_mgof_ee(
 pointer(real scalar function) vector s0, // statistics to be tested
 real colvector f,     // observed count
 real colvector h0,    // expected count
 | real scalar lambda, // lambda for Cressie-Read statistic
   real scalar force,  // do not use partitions (for ksmirnov)
   real scalar dots,   // display dots
   real scalar counter) // will be replaced
{
    real scalar                           n, i, k, nstat, uni, w
    real scalar                           doti, ndots, reps
    real colvector                        fj, j, H, stat, statj, p, hp
    pointer(real colvector) scalar        h
    pointer(real scalar function) vector  s
    struct mm_subsetinfo scalar           info

// set defaults and compute hypothetical distribution
    if (args()<6) dots = 0
    if (s0==NULL) s = &_mm_mgof_mlnp()
    else {
        if (rows(s0)!=1) s = &_mm_mgof_mlnp() \ s0
        else             s = &_mm_mgof_mlnp() , s0
    }
    n = colsum(f)
    k = rows(f)
    uni = (force==1 ? 0 : 1)
    if (rows(h0)==1) {
        h = &(n :* J(k,1,1/k))
        H = rangen(1/k,1,k)
        hp = J(k,1,1/k)
    }
    else {
        h = &h0
        if (mm_isconstant(*h)==0) uni = 0 // reset to 0 if h not constant
        H = mm_colrunsum(*h)
        H = H / H[rows(H)]
        hp = (*h) / colsum(*h)
    }
// compute statistics and p-values
    nstat = length(s)
    stat = statj = p = j = J(nstat,1,0)
    for (i=1; i<=nstat; i++)  stat[i] = (*s[i])(f,*h, lambda, n, hp, H)
    if (uni) {         // uniform null distribution: work with partitions
        if (dots) {
            reps = mm_npartitions(n,k)
            printf("\n{txt}Percent completed ({res}%g{txt} partitions)\n",reps)
            display("{txt}0 {hline 5} 20 {hline 6} 40 {hline 6} 60 {hline 6} 80 {hline 5} 100")
            ndots = (reps>=50 ? 1 : (reps>=25 ? 2 : (reps>=10 ? 5 :
                    (reps>=5 ? 10 : (reps>=2 ? 25 : 50)))))
            doti = 1
        }
        info = _mm_partitionsetup(n, k)
        while ((fj = _mm_partition(info,1)) != J(0,1,.)) {
            w = exp(lnfactorial(info.k) - quadcolsum(lnfactorial(_mm_panels(fj))))
            for (i=1; i<=nstat; i++) {
                statj[i] = (*s[i])(fj,*h, lambda, n, hp, H)
                if (stat[i] <= statj[i]) {
                    j[i] = j[i] + w // statj[1] = -ln(p) => p = exp(-statj[1])
                    p[i] = p[i] + (exp(-statj[1])-p[i])*w/j[i]  // mean updating
                }
            }
            if (dots) {
                if (info.counter*50 >= doti*ndots*reps) {
                    printf(ndots*".")
                    displayflush()
                    doti++
                }
            }
        }
        if (dots) display("")
        counter = info.counter
        return(stat, p:*j)
    }
    if (dots) {
        reps = mm_ncompositions(n,k)
        printf("\n{txt}Percent completed ({res}%g{txt} compositions)\n",reps)
        display("{txt}0 {hline 5} 20 {hline 6} 40 {hline 6} 60 {hline 6} 80 {hline 5} 100")
        ndots = (reps>=50 ? 1 : (reps>=25 ? 2 : (reps>=10 ? 5 :
                (reps>=5 ? 10 : (reps>=2 ? 25 : 50)))))
        doti = 1
    }
    info = _mm_compositionsetup(n, k)
    while ((fj = _mm_composition(info)) != J(0,1,.)) {
        for (i=1; i<=nstat; i++) {
            statj[i] = (*s[i])(fj,*h, lambda, n, hp, H)
            if (stat[i] <= statj[i]) {
                j[i] = j[i] + 1  // statj[1] = -ln(p) => p = exp(-statj[1])
                p[i] = p[i] + (exp(-statj[1])-p[i])/j[i]  // mean updating
            }
        }
        if (dots) {
            if (info.counter*50 >= doti*ndots*reps) {
                printf(ndots*".")
                displayflush()
                doti++
            }
        }
    }
    if (dots) display("")
    counter = info.counter
    return(stat, p:*j)
}

// Pearson's chi2 statistic subroutine
real scalar _mm_mgof_x2(
    real colvector f,           // observed count
    real colvector h,           // expected count
    | transmorphic matrix opt1, // not used
      transmorphic matrix opt2, // not used
      transmorphic matrix opt3, // not used
      transmorphic matrix opt4) // not used
{
    return(quadcolsum((f-h):^2:/h))
}

// log likelihood ratio statistic subroutine
real scalar _mm_mgof_lr(
    real colvector f,           // observed count
    real colvector h,           // expected count
    | transmorphic matrix opt1, // not used
      transmorphic matrix opt2, // not used
      transmorphic matrix opt3, // not used
      transmorphic matrix opt4) // not used
{
    return(2*quadcolsum(f:*(ln(f)-ln(h))))
}

// Cressie-Read statistic
real scalar _mm_mgof_cr(
    real colvector f,           // observed count
    real colvector h,           // expected count
    | real scalar lambda,       // lambda parameter; default is 2/3
      transmorphic matrix opt2, // not used
      transmorphic matrix opt3, // not used
      transmorphic matrix opt4) // not used
{
    real scalar L
    L = (lambda<. ? lambda : 2/3)
    if (abs(L)<1e-6)        return( 2 * quadcolsum(f:*(ln(f)-ln(h))) )
    else if (abs(L+1)<1e-6) return( 2 * quadcolsum(h:*(ln(h)-ln(f))) )
    else if (L>0) return( 2/(L*(L+1)) * quadcolsum(f:*((f:/h):^L :- 1)) )
                  return( 2/(L*(L+1)) * quadcolsum(f:*((h:/f):^-L :- 1)) )
}

// -ln(p) statistic subroutine (multinomial test)
real scalar _mm_mgof_mlnp(
    real colvector f,           // observed count
    real colvector h,           // expected count
    | transmorphic matrix opt1, // not used
      real scalar n,            // n. of obs [ =colsum(f) ]
      real colvector p,         // hypothical propotions [ = h/colsum(h) ]
      transmorphic matrix opt4) // not used
{
    if (args()>=5) return(-(lnfactorial(n) -
     quadcolsum(lnfactorial(f)) + quadcolsum(f:*ln(p))))
    return(-(lnfactorial(quadcolsum(f)) -
     quadcolsum(lnfactorial(f)) + quadcolsum(f:*ln(h/quadcolsum(h)))))
}

// kolmogorov-smirnov statistic subroutine
real scalar _mm_mgof_ksmirnov(
    real colvector f,           // observed count
    real colvector h,           // expected count
    | transmorphic matrix opt1, // not used
      real scalar n,            // n. of obs [ =colsum(f) ]
      transmorphic matrix opt3, // not used
      real colvector H)         // hypothetical cumulative [ = runsum(h)/sum(h) ]
{
    if (args()>=6) return(max(abs(H-mm_colrunsum(f)/n)))
    return(max(abs(mm_colrunsum(h)/colsum(h)-mm_colrunsum(f)/colsum(f))))
}


end
