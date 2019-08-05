*! version 1.0.4  05aug2019  Ben Jann
version 11.2
local Bool  real scalar
local BoolR real rowvector
local Int   real scalar
local IntV  real vector
local IntC  real colvector
local RS    real scalar
local RV    real vector
local RC    real colvector
local RR    real rowvector
local RM    real matrix
local SS    string scalar
local T     transmorphic
local pRC   pointer (`RC') scalar
local pRM   pointer (`RM') scalar
local S     struct _mm_ebal_struct scalar
mata:

struct _mm_ebal_struct {
    `T'    O         // optimization object
    `RS'   N         // target group size (sum of weights)
    `RM'   C         // constraint matrix (excluding collinear columns)
    `RM'   CC        // collinear columns from constraint matrix
    `RC'   Q         // base weights
    `RC'   W         // balancing weights
    `Bool' nc        // normalizing constraint included in optimization
    `RS'   btol      // balancing tolerance
    `Bool' balanced  // 1 if balance achieved, 0 else
    `RR'   Z         // coefficients
    `RR'   g         // gradient
    `RS'   v         // max difference
    `Bool' conv      // 1 if optimize() converged
    `Int'  rc        // return code of optimize()
    `Int'  i         // number of iterations
}

`RR'   mm_ebal_N(`S' S)        return(S.N)
`RM'   mm_ebal_C(`S' S)        return(S.C)
`RM'   mm_ebal_CC(`S' S)       return(S.CC)
`RC'   mm_ebal_Q(`S' S)        return(S.Q)
`RC'   mm_ebal_W(`S' S)        return(S.W)
`Bool' mm_ebal_balanced(`S' S) return(S.balanced)
`RR'   mm_ebal_g(`S' S)        return(S.g)
`RS'   mm_ebal_v(`S' S)        return(S.v)
`Bool' mm_ebal_conv(`S' S)     return(S.conv)
`Int'  mm_ebal_rc(`S' S)       return(S.rc)
`Int'  mm_ebal_i(`S' S)        return(S.i)

`S' mm_ebal_init(
    `RM'   X1,       // treatment group data
    `RC'   w1,       // treatment group base weights
    `RM'   X0,       // control group data
    `RC'   w0,       // control group base weights
  | `IntV' tar,      // balancing targets
    `Bool' cov,      // whether to balance covariances
    `Bool' nc,       // include normalizing constraint in optimization
    `Bool' dfc,      // apply degrees of freedom correction
    `Bool' nostd)    // do not standardize data
{
    `S'     S        // main struct
    `RM'    C1       // expanded X1
    `RM'    M        // target moments
    `RR'    sd       // standard deviations of C1
    `BoolR' collin   // collinear terms
    `RS'    N0       // size of control group (sum of weights)
    
    // defaults
    if (args()<5)  tar   = 1
    if (args()<6)  cov   = 0
    if (args()<7)  nc    = 0
    if (args()<8)  dfc   = 0
    if (args()<9)  nostd = 0
    S.nc   = nc
    S.btol = 1e-5
    // compute group sizes and check input
    // - treatment group
    if (hasmissing(X1)) _error(3351, "missing values in X1 not allowed")
    if (hasmissing(w1)) _error(3351, "missing values in w1 not allowed")
    if (any(w1:<0))     _error(3498, "negative values in w1 not allowed")
    if (rows(w1)==1) S.N = w1*rows(X1)
    else {
        if (rows(w1)!=rows(X1)) _error(3200, "X1 and w1 not conformable")
        S.N = quadsum(w1)
    }
    if (S.N<=0) _error(3498, "treatment group size (sum of weights) must be > 0")
    // - control group
    if (hasmissing(X0)) _error(3351, "missing values in X0 not allowed")
    if (hasmissing(w0)) _error(3351, "missing values in w0 not allowed")
    if (any(w0:<0))     _error(3498, "negative values in w0 not allowed")
    if (rows(w0)==1) N0 = w0*rows(X0)
    else {
        if (rows(w0)!=rows(X0)) _error(3200, "X0 and w0 not conformable")
        N0 = quadsum(w0)
    }
    if (N0<=0) _error(3498, "control group size (sum of weights) must be > 0")
    // - number of variables
    if (cols(X1)!=cols(X0)) _error(3200, "X1 and X0 must contain the same number of columns")
    // expand X1 and X0
    C1  = _mm_ebal_C(tar, cov, X1)
    S.C = _mm_ebal_C(tar, cov, X0)
    // determine collinear columns (across both groups)
    collin = _mm_ebal_collin((C1 \ S.C),
             ((rows(w1)==1 ? J(rows(C1), 1, w1) : w1)
            \ (rows(w0)==1 ? J(rows(S.C), 1, w0) : w0)))
    // compute target moments
    M = mean(C1, w1)
    // apply degrees-of-freedom correction to target moments
    if (dfc) _mm_ebal_dfc(tar, cov, M, cols(X1), (1-1/S.N) / (1-1/N0))
    // remove collinear columns
    if (any(collin)) {
        C1  = select(C1, collin:==0)
        M   = select(M, collin:==0)
        S.C = select(S.C, collin:==0)
    }
    // determine scale for standardization
    if (rows(C1)==1) sd = J(1, cols(C1), 1)
    else {
        if (nostd | cols(C1)==0) sd = J(1, cols(C1), 1)
        else {
            sd = sqrt(diagonal(variance(C1, w1)))'
            _editvalue(sd, 0, 1) // (sd must not be zero)
        }
    }
    // prepare constraints matrix
    // - determine collinear columns
    collin = _mm_ebal_collin(S.C, w0)
    // - standardize constraints
    S.C = (S.C :- M) :/ sd
    // - put collinear columns aside
    if (any(collin)) {
        S.CC = select(S.C, collin)
        S.C  = select(S.C, collin:==0)
    }
    // - add normalization constraint
    if (S.nc) S.C = J(rows(X0),1,1), S.C
    // prepare base weights
    if (S.nc) S.Q = w0
    else      S.Q = w0 / N0
    // prepare optimization object
    S.O = optimize_init()
    optimize_init_which(S.O, "min")
    optimize_init_evaluator(S.O, &_mm_ebal_eval())
    optimize_init_evaluatortype(S.O, "d2")
    optimize_init_conv_ignorenrtol(S.O, "on")
    optimize_init_valueid(S.O, "max difference")
    // set starting values
    if (S.nc) S.Z = (ln(S.N/N0), J(1, cols(S.C)-1, 0))
    else      S.Z = J(1, cols(S.C), 0)
    // done
    return(S)
}

`RM' _mm_ebal_C(`RV' t, `Bool' cov, `RM' X)
{
    `RS' i, j, k, l, c
    `RM' C
    
    // confirm that t only contains 1, 2, or 3
    l = length(t)
    for (i=1; i<=l; i++) {
        if (!anyof((1,2,3), t[i])) _error("invalid target")
    }
    // means
    k = cols(X)
    C = X
    // variances
    if (any(t:>=2)) {
        c = 0
        for (i=1; i<=k; i++) { // count number of terms
            if (t[mod(i-1, l) + 1] >= 2) c++
        }
        C = C, J(rows(C), c, .)
        c = cols(C) - c
        for (i=1; i<=k; i++) {
            if (t[mod(i-1, l) + 1] >= 2) {
                C[,++c] = X[,i]:^2
            }
        }
    }
    // covariances
    if (cov) {
        c = (k^2 - k)/2
        C = C, J(rows(C), c, .)
        c = cols(C) - c
        for (i=1; i<k; i++) {
            for (j=i+1; j<=k; j++) {
                C[,++c] = X[,i] :* X[,j]
            }
        }
    }
    // skewnesses
    if (anyof(t,3)) {
        c = 0
        for (i=1; i<=k; i++) { // count number of terms
            if (t[mod(i-1, l) + 1] == 3) c++
        }
        C = C, J(rows(C), c, .)
        c = cols(C) - c
        for (i=1; i<=k; i++) {
            if (t[mod(i-1, l) + 1] == 3) {
                C[,++c] = X[,i]:^3
            }
        }
    }
    return(C)
}

`BoolR' _mm_ebal_collin(`RM' C, `RC' w)
{
    `RR' m
    `RM' CP

    if (cols(C)==0) return(J(1,0,.))
    m  = mean(C, w)
    CP = quadcrossdev(C, m, w, C, m)
    return((diagonal(invsym(CP, 1..cols(CP))):==0)')
}

void _mm_ebal_dfc(`RV' t, `Bool' cov, `RM' M, `RS' k, `RS' dfc)
{
    `RS' i, j, l, c
    `RM' V
    
    if (cols(M)==k) return
    if (dfc==1)     return
    c = k
    l = length(t)
    // variances
    if (any(t:>=2)) {
        V = J(2, k, 1) // needed for skewness correction
        for (i=1; i<=k; i++) {
            if (t[mod(i-1, l) + 1] >= 2) {
                c++
                V[1,i] = M[c]
                M[c]   = (M[c] - (1-dfc) * M[i]^2) / dfc
                V[2,i] = M[c]
            }
        }
    }
    // covariances
    if (cov) {
        for (i=1; i<k; i++) {
            for (j=i+1; j<=k; j++) {
                c++
                M[c] = (M[c] - (1-dfc) * M[i] * M[j]) / dfc
            }
        }
    }
    // skewnesses
    if (anyof(t,3)) {
        for (i=1; i<=k; i++) {
            if (t[mod(i-1, l) + 1] == 3) {
                c++
                M[c] = (M[c] - (V[1,i] - dfc^1.5 * V[2,i]) * 3 *M[i]
                             + (1 - dfc^1.5) * 2 * M[i]^3) / dfc^1.5
            }
        }
    }
}

`T' mm_ebal_btol(`S' S, | `RS' btol)
{
    if (args()==1) return(S.btol)
    S.btol = btol
}

`T' mm_ebal_trace(`S' S, | `SS' tracelevel)
{
    if (args()==1) return(optimize_init_tracelevel(S.O))
    optimize_init_tracelevel(S.O, tracelevel)
}

`T' mm_ebal_difficult(`S' S, | `Bool' flag)
{
    if (args()==1) return(optimize_init_singularHmethod(S.O)=="hybrid")
    if (flag)      optimize_init_singularHmethod(S.O, "hybrid")
    else           optimize_init_singularHmethod(S.O, "")
}

`T' mm_ebal_maxiter(`S' S, | `Int' maxiter)
{
    if (args()==1) return(optimize_init_conv_maxiter(S.O))
    optimize_init_conv_maxiter(S.O, maxiter)
}

`T' mm_ebal_ptol(`S' S, | `RS' tol)
{
    if (args()==1) return(optimize_init_conv_ptol(S.O))
    optimize_init_conv_ptol(S.O, tol)
}

`T' mm_ebal_vtol(`S' S, | `RS' tol)
{
    if (args()==1) return(optimize_init_conv_vtol(S.O))
    optimize_init_conv_vtol(S.O, tol)
}

`T' mm_ebal_nowarn(`S' S, | `Bool' flag)
{
    if (args()==1) return(optimize_init_conv_warning(S.O)=="off")
    if (flag)      optimize_init_conv_warning(S.O, "off")
    else           optimize_init_conv_warning(S.O, "on")
}

`T' mm_ebal_Z(`S' S, | `RR' Z)
{
    if (args()==1) return(S.Z)
    S.Z = Z
}

`Bool' mm_ebal(`S' S)
{
    // find solution
    if (S.nc==0 & cols(S.C)==0) {
        // no covariates; nothing to do
        S.Z = S.g = J(1,0,.)
        S.conv = 1
        S.rc = S.i = S.v = 0
    }
    else {
        optimize_init_argument(S.O, 1, S.nc)
        optimize_init_argument(S.O, 2, &S.Q)
        optimize_init_argument(S.O, 3, &S.C)
        optimize_init_argument(S.O, 4, S.N)
        optimize_init_params(S.O, S.Z)
        (void) _optimize(S.O)
        if (optimize_result_errorcode(S.O)) {
            errprintf("{p}\n")
            errprintf("%s\n", optimize_result_errortext(S.O))
            errprintf("{p_end}\n")
        }
        S.Z    = optimize_result_params(S.O)
        S.conv = optimize_result_converged(S.O)
        S.rc   = optimize_result_returncode(S.O)
        S.i    = optimize_result_iterations(S.O)
        S.v    = optimize_result_value(S.O)
        S.g    = optimize_result_gradient(S.O)
    }
    // obtain balancing weights from solution and update balancing loss
    if (S.nc) {
        S.W = S.Q :* exp(quadcross(S.C', S.Z'))
        if (cols(S.CC)>0) S.v = max((S.v, abs(quadcross(S.W, S.CC)) / S.N))
    }
    else {
        S.W = quadcross(S.C', S.Z')
        S.W = S.Q :* exp(S.W:-max(S.W))
        S.W = S.W / quadsum(S.W)
        if (cols(S.CC)>0) S.v = max((S.v, abs(quadcross(S.W, S.CC))))
        S.W = S.W * S.N
    }
    if (cols(S.CC)>0) {
        if (optimize_init_tracelevel(S.O)!="none") {
            printf("{txt}Final fit including collinear terms:\n")
            printf("               max difference = {res}%10.0g\n", S.v)
        }
    }
    // check whether balancing criterion fulfilled
    S.balanced = (S.v < S.btol)
    return(S.balanced)
}

void _mm_ebal_eval(`Int' todo, `RR' Z, 
    `Bool' nc, `pRC' Q, `pRM' C, `RS' N, 
    `RS' v, `RR' g, `RM' H)
{
    `RC' W
    
    if (nc) {
        W = *Q :* exp(quadcross(*C', Z'))
        g = quadcross(W, *C)
        g[1] = g[1] - N
        v = max(abs(g)) / N
    }
    else {
        W = quadcross(*C', Z')
        W = *Q :* exp(W:-max(W)) // set exp(max)=1 to avoid numerical overflow
        W = W / quadsum(W)
        g = quadcross(W, *C)
        v = max(abs(g))
    }
    if (todo==2) {
        H = quadcross(*C, W, *C)
    }
}

end

