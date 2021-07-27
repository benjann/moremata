*! version 1.0.0  27jul2021  Ben Jann

version 11.2

// class & struct
local MAIN   mm_ebalance
local SETUP  _`MAIN'_setup
local Setup  struct `SETUP' scalar
local IF     _`MAIN'_IF
local If     struct `IF' scalar
// real
local RS     real scalar
local RR     real rowvector
local RC     real colvector
local RM     real matrix
// counters
local Int    real scalar
local IntC   real colvector
// string
local SS     string scalar
// boolean
local Bool   real scalar
local BoolC  real colvector
// transmorphic
local T      transmorphic
// pointers
local PC     pointer(real colvector) scalar
local PM     pointer(real matrix) scalar

mata:

// class ----------------------------------------------------------------------

struct `SETUP' {
    // data
    `PM'    X, X0    // pointers to main data and reference data
    `PC'    w, w0    // pointers to base weights
    `Int'   N, N0    // number of obs
    `RS'    W, W0    // sum of weights
    `RR'    m, m0    // means
    `Int'   k        // number of terms
    `BoolC' omit     // flag collinear terms
    `Int'   k_omit   // number of omitted terms
    
    // settings
    `SS'   ltype      // type of loss function
    `SS'   trace      // trace level
    `Bool' difficult  // use hybrid optimization
    `Int'  maxiter    // max number of iterations
    `RS'   ptol       // convergence tolerance for the parameter vector
    `RS'   vtol       // convergence tolerance for the balancing loss
    `RC'   btol       // balancing tolerance
    `Bool' nowarn     // do not display no convergence/balance warning
}

struct `IF' {
    `RM'   b, b0      // influence functions of coefficients
    `RC'   a, a0      // influence function of intercept
}

class `MAIN' {
    // settings
    private:
        void    new()           // initialize class with default settings
        void    clear()         // clear all results
        `Bool'  nodata()        // whether data is set
        `Setup' setup           // container for data and settings
    public:
        void    data()
        `RM'    X(), Xref()     // retrieve data
        `RC'    w(), wref()     // retrieve base weights
        `RS'    N(), Nref()     // retrieve number of obs
        `RS'    W(), Wref()     // retrieve sum of weights
        `RR'    m(), mref()     // retrieve moments
        `Int'   k()             // retrieve number of terms
        `RC'    omit()          // retrieve omitted flags
        `Int'   k_omit()        // retrieve number of omitted terms
        `T'     ltype()         // set/retrieve loss function
        `T'     trace()         // set/retrieve trace level
        `T'     difficult()     // set/retrieve difficult flag
        `T'     maxiter()       // set/retrieve max iterations
        `T'     ptol()          // set/retrieve p-tolerance
        `T'     vtol()          // set/retrieve v-tolerance
        `T'     btol()          // set/retrieve balancing tolerance
        `T'     nowarn()        // set/retrieve nowarn flag
    
    // results
    public:
        `RC'    b()             // retrieve coefficients
        `RS'    a()             // retrieve normalizing intercept
        `RC'    wbal()          // retrieve balancing weights
        `RR'    madj()          // adjusted (reweighted) means
        `Int'   iter()          // retrieve number of iterations
        `Bool'  converged()     // retrieve convergence flag
        `RS'    loss()          // retrieve final balancing loss
        `Bool'  balanced()      // retrieve balancing flag
        `RM'    IF_b(), IFref_b() // retrieve IF of coefficients
        `RC'    IF_a(), IFref_a() // retrieve IF of intercept
    private:
        `RC'    b               // coefficients
        `RS'    a               // normalizing intercept
        `RC'    wbal            // balancing weights
        `RR'    madj            // adjusted (reweighted) means
        `Int'   iter            // number of iterations
        `Bool'  conv            // optimize() convergence
        `RC'    loss            // balancing loss
        `Bool'  balanced        // balance achieved
        `If'    IF              // influence functions
        void    Fit(), _Fit_a(), _Fit_b() // fit coefficients
        void    _IF_b(), _IF_a() // generate influence functions
}

// init -----------------------------------------------------------------------

void `MAIN'::new()
{
    setup.ltype     = "reldif"
    setup.trace     = (st_global("c(iterlog)")=="off" ? "none" : "value")
    setup.difficult = 0
    setup.maxiter   = st_numscalar("c(maxiter)")
    setup.ptol      = 1e-6
    setup.vtol      = 1e-7
    setup.btol      = 1e-6
    setup.nowarn    = 0
}

void `MAIN'::clear()
{
    b    = J(0,1,.)
    a    = .
    wbal = J(0,1,.)
    madj = J(1,0,.)
    iter = conv = loss = balanced = .
    IF   = `IF'()
}

// data -----------------------------------------------------------------------

void `MAIN'::data(`RM' X, `RC' w, `RM' X0, `RC' w0, | `Bool' fast)
{
    `RM' CP
    
    // check for missing values and negative weights
    if (args()<5) fast = 0
    if (!fast) {
        if (missing(X) | missing(w) | missing(X0) | missing(w0)) _error(3351)
        if (any(w:<0) | any(w0:<0)) _error(3498, "w and wref must be positive")
    }
    // obtain main data
    setup.X = &X
    setup.w = &w
    setup.N = rows(X)
    if (setup.N==0) _error(2000, "no observations in main data")
    if (rows(w)!=1) {
        if (rows(X)!=rows(w)) _error(3200, "X and w not conformable")
        setup.W = quadsum(w)
    }
    else setup.W = setup.N * w
    if (setup.N!=0 & setup.W==0) _error(3498, "sum(w) must be > 0")
    // obtain reference data
    setup.X0 = &X0
    setup.w0 = &w0
    setup.N0 = rows(X0)
    if (setup.N0==0) _error(2000, "no d observations in reference data")
    if (rows(w0)!=1) {
        if (rows(X0)!=rows(w0)) _error(3200, "Xref and wref not conformable")
        setup.W0 = quadsum(w0)
    }
    else setup.W0 = setup.N0 * w0
    if (setup.N0!=0 & setup.W0==0) _error(3498, "sum(wref) must be > 0")
    // target moments
    setup.k = cols(X)
    if (setup.k!=cols(X0)) _error(3200, "X and Xref not conformable")
    setup.m0 = mean(X0, w0)
    // identify collinear terms
    setup.m  = mean(X, w)
    if (setup.k==0) {
        setup.omit = J(0,1,.)
        setup.k_omit = 0
    }
    else {
        CP = quadcrossdev(X, setup.m, w, X, setup.m)
        setup.omit = (diagonal(invsym(CP)):==0) // or: diagonal(invsym(CP, 1..setup.k)):==0
        setup.k_omit = sum(setup.omit)
    }
    // clear results
    clear()
}

`Bool' `MAIN'::nodata() return(setup.X==NULL)

`RM' `MAIN'::X()
{
    if (nodata()) return(J(0,1,.))
    return(*setup.X)
}

`RM' `MAIN'::Xref()
{
    if (nodata()) return(J(0,1,.))
    return(*setup.X0)
}

`RC' `MAIN'::w()
{
    if (nodata()) return(J(0,1,.))
    return(*setup.w)
}

`RC' `MAIN'::wref()
{
    if (nodata()) return(J(0,1,.))
    return(*setup.w0)
}

`Int' `MAIN'::N() return(setup.N)

`Int' `MAIN'::Nref() return(setup.N0)

`RS' `MAIN'::W() return(setup.W)

`RS' `MAIN'::Wref() return(setup.W0)

`RR' `MAIN'::m() return(setup.m)

`RR' `MAIN'::mref() return(setup.m0)

`Int' `MAIN'::k() return(setup.k)

`RC'  `MAIN'::omit() return(setup.omit)

`Int' `MAIN'::k_omit() return(setup.k_omit)

// settings -------------------------------------------------------------------

`T' `MAIN'::ltype(| `SS' ltype)
{
    if (args()==0) return(setup.ltype)
    if (setup.ltype==ltype) return // no change
    if (!anyof(("reldif", "absdif"), ltype)) {
        printf("{err}'%s' not allowed\n", ltype)
        _error(3498)
    }
    setup.ltype = ltype
    clear()
}

`T' `MAIN'::trace(| `SS' trace)
{
    `T' S
    
    if (args()==0) return(setup.trace)
    if (setup.trace==trace) return // no change
    S = optimize_init()
    optimize_init_tracelevel(S, trace) // throw error if trace is invalid
    setup.trace = trace
    clear()
}

`T' `MAIN'::difficult(| `Bool' difficult)
{
    if (args()==0) return(setup.difficult)
    if (setup.difficult==(difficult!=0)) return // no change
    setup.difficult = (difficult!=0)
    clear()
}

`T' `MAIN'::maxiter(| `Int' maxiter)
{
    if (args()==0) return(setup.maxiter)
    if (setup.maxiter==maxiter) return // no change
    if (maxiter<0) _error(3498, "setting out of range")
    setup.maxiter = maxiter
    clear()
}

`T' `MAIN'::ptol(| `RS' ptol)
{
    if (args()==0) return(setup.ptol)
    if (setup.ptol==ptol) return // no change
    if (ptol<=0) _error(3498, "setting out of range")
    setup.ptol = ptol
    clear()
}

`T' `MAIN'::vtol(| `RS' vtol)
{
    if (args()==0) return(setup.vtol)
    if (setup.vtol==vtol) return // no change
    if (vtol<=0) _error(3498, "setting out of range")
    setup.vtol = vtol
    clear()
}

`T' `MAIN'::btol(| `RS' btol)
{
    if (args()==0) return(setup.btol)
    if (setup.btol==btol) return // no change
    if (btol<=0) _error(3498, "setting out of range")
    setup.btol = btol
    clear()
}

`T' `MAIN'::nowarn(| `Bool' nowarn)
{
    if (args()==0) return(setup.nowarn)
    if (setup.nowarn==(nowarn!=0)) return // no change
    setup.nowarn = (nowarn!=0)
    clear()
}

// results --------------------------------------------------------------------

`RC' `MAIN'::b()
{
    if (length(b)==0) Fit()
    return(b)
}

`RS' `MAIN'::a()
{
    if (length(b)==0) Fit()
    return(a)
}

`RC' `MAIN'::wbal()
{
    if (length(b)==0) Fit()
    return(wbal)
}

`RR' `MAIN'::madj()
{
    if (length(madj)) return(madj)
    if (length(b)==0) Fit()
    madj = mean(X(), wbal)
    return(madj)
}

`Int' `MAIN'::iter()
{
    
    if (length(b)==0) Fit()
    return(iter)
}

`Bool' `MAIN'::converged()
{
    if (length(b)==0) Fit()
    return(conv)
}

`RS' `MAIN'::loss()
{
    if (length(b)==0) Fit()
    return(loss)
}

`RS' `MAIN'::balanced()
{
    if (length(b)==0) Fit()
    return(balanced)
}

`RM' `MAIN'::IF_b()
{
    if (rows(IF.b)==0) _IF_b()
    return(IF.b)
}

`RM' `MAIN'::IFref_b()
{
    if (rows(IF.b0)==0) _IF_b()
    return(IF.b0)
}

`RC' `MAIN'::IF_a()
{
    if (rows(IF.a)==0) _IF_a()
    return(IF.a)
}

`RC' `MAIN'::IFref_a()
{
    if (rows(IF.a0)==0) _IF_a()
    return(IF.a0)
}

// optimization ---------------------------------------------------------------

void `MAIN'::Fit()
{
    // optimize
    if (nodata()) _error(3498, "data not set")
    if ((k()-k_omit())<=0) { // no covariates
        b = J(k(), 1, 0)
        iter = loss = 0
        conv = 1
    }
    else _Fit_b() // fit coefficients
    _Fit_a()      // compute intercept (and balancing weights)
    
    // check balancing
    if (k_omit()) {
        // recompute balancing loss in case of omitted terms
        loss = _mm_ebalance_loss(ltype(), mean(X():-mref(), wbal), mref())
    }
    if (trace()!="none") {
        printf("{txt}Final fit:     balancing loss = {res}%10.0g\n", loss)
    }
    balanced = (loss<btol())
    if (balanced) return
    if (nowarn()) return
    display("{txt}balance not achieved")
}

void `MAIN'::_Fit_a()
{
    `RS' ul
    `RC' xb
    
    xb = X() * b
    ul = max(xb) // set exp(max)=1 to avoid numerical overflow
    a  = ln(Wref()) - ln(quadsum(w() :* exp(xb :- ul))) - ul
    wbal = w() :* exp(xb :+ a)
}

void `MAIN'::_Fit_b()
{
    `IntC' p
    `T'    S
    
    // setup
    if (k_omit()) p = select(1::k(), omit():==0)
    S = optimize_init()
    optimize_init_which(S, "min")
    optimize_init_evaluator(S, &_mm_ebalance_eval())
    optimize_init_evaluatortype(S, "d2")
    optimize_init_technique(S, "nr")
    optimize_init_singularHmethod(S, difficult() ? "hybrid" : "")
    optimize_init_conv_maxiter(S, maxiter())
    optimize_init_conv_ptol(S, ptol())
    optimize_init_conv_vtol(S, vtol())
    optimize_init_conv_ignorenrtol(S, "on")
    optimize_init_conv_warning(S, nowarn() ? "off" : "on")
    optimize_init_tracelevel(S, trace())
    optimize_init_valueid(S, "balancing loss")
    optimize_init_params(S, J(1, k()-k_omit(), 0))
    optimize_init_argument(S, 1, k_omit() ? X()[,p]:-mref()[p] : X():-mref())
    optimize_init_argument(S, 2, w())
    optimize_init_argument(S, 3, k_omit() ? mref()[p] : mref())
    optimize_init_argument(S, 4, ltype())
    
    // run optimizer
    (void) _optimize(S)
    if (optimize_result_errorcode(S)) {
        errprintf("{p}\n")
        errprintf("%s\n", optimize_result_errortext(S))
        errprintf("{p_end}\n")
        exit(optimize_result_returncode(S))
    }
    
    // obtain results
    if (k_omit()) {
        b = J(k(), 1, 0)
        b[p] = optimize_result_params(S)'
    }
    else b = optimize_result_params(S)'
    iter = optimize_result_iterations(S)
    loss = optimize_result_value(S)
    conv = optimize_result_converged(S)
}

void _mm_ebalance_eval(`Int' todo, `RR' b, `RM' X, `RC' w0, `RR' m, `SS' ltype,
    `RS' v, `RR' g, `RM' H)
{
    `RS' W
    `RC' w
    
    w = X * b'
    w = w0 :* exp(w :- max(w)) // set exp(max)=1 to avoid numerical overflow
    W = quadsum(w)
    g = quadcross(w, X) / W
    v = _mm_ebalance_loss(ltype, g, m) //mreldif(g+m, m)
    if (todo==2) H = quadcross(X, w, X) / W
}

`RS' _mm_ebalance_loss(`SS' ltype, `RR' d, `RR' m)
{
    if (ltype=="absdif") return(max(abs(d)))
    return(mreldif(d+m, m)) // ltype=="reldif"
}

// influence functions --------------------------------------------------------

void `MAIN'::_IF_b()
{
    if (length(b)==0) Fit()
    _mm_ebalance_IF_b(IF, X(), Xref(), w(), wbal, mref())
}

void `MAIN'::_IF_a()
{
    if (rows(IF.b)==0) _IF_b()
    _mm_ebalance_IF_a(IF, madj(), w(), wbal, Wref())
}

void _mm_ebalance_IF_b(`If' IF, `RM' X, `RM' Xref, `RC' w, `RC' wbal, `RR' mref)
{
    `RM' Q
    
    IF.b = X :- mref                        // moment condition of b
    Q = -invsym(quadcross(IF.b, wbal, X))   // derivative
    IF.b  = wbal:/w :* (IF.b * Q')          // contribution of main sample
    IF.b0 = (mref :- Xref) * Q'             // contribution of reference sample
}

void _mm_ebalance_IF_a(`If' IF, `RR' madj, `RC' w, `RC' wbal, `RS' Wref)
{
    `RM' Q
    
    Q = madj * Wref                         // cross-derivative
    IF.a  = -(wbal:/w :+ IF.b * Q') / Wref  // contribution of main sample
    IF.a0 =  (1 :- IF.b0 * Q') / Wref       // contribution of reference sample
}

end

