*! version 1.0.1  22mar2021  Ben Jann
version 9.2
mata:

real colvector mm_lsfit(real colvector y, | real matrix X, real colvector w,
    real scalar cons, real scalar qd, real scalar demean)
{
    if (args()<3) w = 1
    return(mm_ls_b(mm_ls(y, X, w, cons, qd, demean)))
}

struct mm_ls_struct {
    pointer(real colvector) scalar y, w
    pointer(real matrix) scalar    X
    real scalar    cons
    real scalar    k    // number of predictors
    real scalar    kadj // number of non-collinear predictors
    real scalar    N    // sum of weights
    real scalar    ymean
    real rowvector means
    real scalar    rss // residual sum of squares
    real scalar    s   // scale (RMSE)
    real scalar    r2  // R-squared
    real matrix    SSinv, XXinv, V
    real colvector omit
    real colvector b
}

struct mm_ls_struct scalar mm_ls(
    real colvector   y,
    | real matrix    X,
      real colvector w,
      real scalar    cons,   // cons=0 excludes the constant
      real scalar    qd,     // qd=0 uses single precision for cross products
      real scalar    demean) // demean=0 does not use demeaning
{
    real scalar     ymean
    real rowvector  means
    real colvector  ybar, p, Xy, beta
    real matrix     Xbar, XX
    struct mm_ls_struct scalar t
    
    // initialize
    if (args()<3) w = 1
    t.y = &y
    t.X = &X
    t.w = &w
    t.cons = (cons!=0)
    t.k = cols(X)
    t.N = t.ymean = t.means = t.rss = t.s = t.r2 = t.XXinv = t.V = .z
    if (!cons & !t.k) { // model has no parameters
        t.b = t.omit = J(0,1,.)
        t.kadj = 0
        t.XXinv = J(0,0,.)
        t.rss = .
        return(t)
    }
    if (rows(y)==0) {
        t.b = t.omit = J(t.k+t.cons,1,.)
        t.XXinv = J(t.k+t.cons,t.k+t.cons,.)
        t.rss = .
        return(t)
    }
    
    // without constant
    if (cons==0) {
        XX = qd ? quadcross(X, w, X) : cross(X, w, X)
        t.XXinv = invsym(XX)
        t.omit = diagonal(t.XXinv):==0
        p = select(1::t.k, !t.omit)
        t.kadj = length(p)
        if (t.kadj==0) {
            t.b = J(t.k,1,0)
            return(t)
        }
        if (t.kadj<t.k) {
            XX = XX[p,p]
            Xy = qd ? quadcross(X[,p], w, y) : cross(X[,p], w, y)
        }
        else Xy = qd ? quadcross(X, w, y) : cross(X, w, y)
        beta = lusolve(XX, Xy)
        t.b = J(t.k, 1, 0)
        t.b[p] = beta
        return(t)
    }
    
    // constant only
    if (t.k==0) {
        t.kadj = t.omit = 0
        t.b = mm_ls_ymean(t)
        return(t)
    }
    
    // without demeaning
    if (!demean) {
        XX = qd ? quadcross(X,1, w, X,1) : cross(X,1, w, X,1)
        t.XXinv = invsym(XX, t.k+1) // do not omit constant
        t.omit = diagonal(t.XXinv)[|1\t.k|]:==0
        p = select(1::t.k, !t.omit)
        t.kadj = length(p)
        if (t.kadj==0) { // all collinear
            t.b = J(t.k,1,0) \ mm_ls_ymean(t) // add constant
            t.omit = t.omit \ 0
            return(t)
        }
        if (t.kadj<t.k) {
            XX = XX[p \ t.k+1, p \ t.k+1]
            Xy = qd ? quadcross(X[,p],1, w, y,0) : cross(X[,p],1, w, y,0)
        }
        else Xy = qd ? quadcross(X,1, w, y,0) : cross(X,1, w, y,0)
        beta = lusolve(XX, Xy)
        t.b = J(t.k+1, 1, 0)
        t.b[p \ t.k+1] = beta
        t.omit = t.omit \ 0 // add constant
        return(t)
    }
    
    // mean deviation method
    ymean = mm_ls_ymean(t)
    means = mm_ls_means(t)
    Xbar = X :- means
    XX = qd ? quadcross(Xbar, w, Xbar) : cross(Xbar, w, Xbar)
    t.SSinv = invsym(XX)
    t.omit = diagonal(t.SSinv):==0
    p = select(1::t.k, !t.omit)
    t.kadj = length(p)
    if (t.kadj==0) { // all collinear
        t.b = J(t.k,1,0) \ ymean // add constant
        t.omit = t.omit \ 0
        return(t)
    }
    if (t.kadj<t.k) {
        Xbar  = Xbar[,p]
        XX    = XX[p,p]
        means = means[p]
    }
    ybar = y :- ymean
    Xy = qd ? quadcross(Xbar, w, ybar) : cross(Xbar, w, ybar)
    beta = lusolve(XX, Xy)
    t.b = J(t.k+1, 1, 0)
    t.b[p]   = beta
    t.b[t.k+1] = (ymean - means*beta) // add constant
    t.omit = t.omit \ 0
    return(t)
}

real colvector mm_ls_b(struct mm_ls_struct scalar t)
{
    return(t.b)
}

real colvector mm_ls_omit(struct mm_ls_struct scalar t)
{
    return(t.omit)
}

real colvector mm_ls_k_omit(struct mm_ls_struct scalar t)
{
    return(t.k - t.kadj)
}

real matrix mm_ls_V(struct mm_ls_struct scalar t)
{
    if (t.V==.z) t.V = mm_ls_XXinv(t)*mm_ls_s(t)^2
    return(t.V)
}

real colvector mm_ls_se(struct mm_ls_struct scalar t)
{
    return(sqrt(diagonal(mm_ls_V(t))))
}

real scalar mm_ls_N(struct mm_ls_struct scalar t)
{
    if (t.N==.z) t.N = rows(*t.w)==1 ? *t.w * rows(*t.y) : quadsum(*t.w)
    return(t.N)
}

real scalar mm_ls_ymean(struct mm_ls_struct scalar t)
{
    if (t.ymean==.z) {
        if (rows(*t.y)==0) t.ymean = .
        else               t.ymean = quadcross(*t.w, *t.y) / mm_ls_N(t)
    }
    return(t.ymean)
}

real rowvector mm_ls_means(struct mm_ls_struct scalar t)
{
    if (t.means==.z) {
        if      (t.k==0)        t.means = J(1, 0, .)
        else if (rows(*t.X)==0) t.means = J(1, t.k, .)
        else                    t.means = quadcross(*t.w, *t.X) / mm_ls_N(t)
    }
    return(t.means)
}

real matrix mm_ls_XXinv(struct mm_ls_struct scalar t)
{
    // note: XXinv will already be filled in case of cons=0 or demean=0
    if (t.XXinv==.z) {
        if (t.k) {
            t.XXinv = t.SSinv \ -(t.means * t.SSinv)
            t.XXinv = (t.XXinv, (t.XXinv[t.k+1,]' \
                       1/mm_ls_N(t) - t.means * t.XXinv[t.k+1,]'))
        }
        else t.XXinv = 1/mm_ls_N(t) // constant-only model
    }
    return(t.XXinv)
}

real colvector mm_ls_xb(struct mm_ls_struct scalar t, | real matrix X)
{
    if (args()==2) {
        if (cols(X)!=t.k) _error(3200)
        if (t.k==0) {
            if (t.cons) return(J(rows(X), 1, t.b))
            return(J(rows(X), 1, .)) // model has no parameters
        }
        if (t.cons) return(X * t.b[|1\t.k|] :+ t.b[t.k+1])
        return(X * t.b)
    }
    if (t.k==0) {
        if (t.cons) return(J(rows(*t.y), 1, t.b))
        return(J(rows(*t.y), 1, .)) // model has no parameters
    }
    if (t.cons) return(*t.X * t.b[|1\t.k|] :+ t.b[t.k+1])
    return(*t.X * t.b)
}

real scalar mm_ls_rss(struct mm_ls_struct scalar t)
{
    if (t.rss==.z) {
        t.rss = quadcross(*t.w, (*t.y - mm_ls_xb(t)):^2)
    }
    return(t.rss)
}

real scalar mm_ls_s(struct mm_ls_struct scalar t)
{
    if (t.s==.z) {
        t.s = sqrt(mm_ls_rss(t) * editmissing(1 / (mm_ls_N(t)-t.kadj-t.cons), 0))
    }
    return(t.s)
}

real scalar mm_ls_r2(struct mm_ls_struct scalar t)
{
    if (t.r2==.z) {
        t.r2 = 1 - mm_ls_rss(t)/quadcross(*t.w, (*t.y :- mm_ls_ymean(t)):^2)
    }
    return(t.r2)
}

end
