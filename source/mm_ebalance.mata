*! version 1.0.8  27apr2022  Ben Jann

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
local IntR   real rowvector
local IntM   real matrix
// string
local SS     string scalar
// boolean
local Bool   real scalar
local BoolC  real colvector
// transmorphic
local T      transmorphic
local TS     transmorphic scalar
// pointers
local PC     pointer(real colvector) scalar
local PM     pointer(real matrix) scalar

mata:

// class ----------------------------------------------------------------------

struct `SETUP' {
    // data
    `PM'    X, X0      // pointers to main data and reference data
    `PC'    w, w0      // pointers to base weights
    `Int'   N, N0      // number of obs
    `RS'    W, W0      // sum of weights
    `RR'    m, m0      // means
    `RR'    s, s0      // scales
    `Int'   k          // number of terms
    `BoolC' omit       // flag collinear terms
    `Int'   k_omit     // number of omitted terms
    
    // settings
    `IntR'  adj, noadj // indices of columns to be adjusted/not adjusted
    `T'     tau        // target sum of weights
    `SS'    scale      // type of scales
    `RC'    btol       // balancing tolerance
    `SS'    ltype      // type of loss function
    `SS'    etype      // evaluator type
    `SS'    trace      // trace level
    `Int'   maxiter    // max number of iterations
    `RS'    ptol       // convergence tolerance for the parameter vector
    `RS'    vtol       // convergence tolerance for the balancing loss
    `Bool'  difficult  // use hybrid optimization
    `Bool'  nostd      // do not standardize
    `Bool'  nowarn     // do not display no convergence/balance warning
}

struct `IF' {
    `RM'   b, b0      // influence functions of coefficients
    `RC'   a, a0      // influence function of intercept
}

class `MAIN' {
    // settings
    private:
        void    new()             // initialize class with default settings
        void    clear()           // clear all results
        `Bool'  nodata()          // whether data is set
        `Setup' setup             // container for data and settings
    public:
        void    data()
        `RM'    X(), Xref()       // retrieve data
        `RC'    w(), wref()       // retrieve base weights
        `RS'    N(), Nref()       // retrieve number of obs
        `RS'    W(), Wref()       // retrieve sum of weights
        `RR'    m(), mref()       // retrieve moments
        `RR'    s(), sref()       // retrieve scales
        `RR'    mu()              // retrieve target moments
        `Int'   k()               // retrieve number of terms
        `RC'    omit()            // retrieve omitted flags
        `Int'   k_omit()          // retrieve number of omitted terms
        `T'     adj(), noadj()    // set/retrieve adj/noadj
        `T'     tau()             // set/retrieve target sum of weights
        `T'     scale()           // set/retrieve scales
        `T'     btol()            // set/retrieve balancing tolerance
        `T'     ltype()           // set/retrieve loss function
        `T'     etype()           // set/retrieve evaluator type
        `T'     alteval()         // set/retrieve alteval flag (old)
        `T'     trace()           // set/retrieve trace level
        `T'     maxiter()         // set/retrieve max iterations
        `T'     ptol()            // set/retrieve p-tolerance
        `T'     vtol()            // set/retrieve v-tolerance
        `T'     difficult()       // set/retrieve difficult flag
        `T'     nostd()           // set/retrieve nostd flag
        `T'     nowarn()          // set/retrieve nowarn flag
    
    // results
    public:
        `RC'    b()               // retrieve coefficients
        `RS'    a()               // retrieve normalizing intercept
        `RC'    xb()              // retrieve linear prediction
        `RC'    wbal()            // retrieve balancing weights
        `RC'    pr()              // retrieve propensity score
        `RR'    madj()            // adjusted (reweighted) means
        `RS'    wsum()            // retrieve sum of balancing weights
        `RS'    loss()            // retrieve final balancing loss
        `Bool'  balanced()        // retrieve balancing flag
        `RS'    value()           // retrieve value of optimization criterion
        `Int'   iter()            // retrieve number of iterations
        `Bool'  converged()       // retrieve convergence flag
        `RM'    IF_b(), IFref_b() // retrieve IF of coefficients
        `RC'    IF_a(), IFref_a() // retrieve IF of intercept
    private:
        `RS'    tau               // target sum of weight
        `RR'    mu                // target means
        `IntM'  adj, noadj        // permutation vectors for source of mu
        `RR'    scale             // scales for standardization
        `RC'    b                 // coefficients
        `RS'    a                 // normalizing intercept
        `RC'    xb                // linear prediction (without a)
        `RC'    wbal              // balancing weights
        `RR'    madj              // adjusted (reweighted) means
        `RS'    wsum              // sum of balancing weights
        `RC'    loss              // balancing loss
        `Bool'  balanced          // balance achieved
        `RS'    value             // value of optimization criterion
        `Int'   iter              // number of iterations
        `Bool'  conv              // optimize() convergence
        `If'    IF                // influence functions
        void    _IF_b(), _IF_a()  // generate influence functions
        void    Fit()             // fit coefficients
        void    _Fit_b(), _Fit_a()
        `RM'    _Fit_b_X(), _Fit_b_Xc()
        `RR'    _Fit_b_mu()
        void    _setadj()         // fill in adj and noadj
}

// init -----------------------------------------------------------------------

void `MAIN'::new()
{
    setup.tau       = "Wref"
    setup.scale     = "main"
    setup.ltype     = "reldif"
    setup.etype     = "bl"
    setup.nostd     = 0
    setup.trace     = (st_global("c(iterlog)")=="off" ? "none" : "value")
    setup.difficult = 0
    setup.maxiter   = st_numscalar("c(maxiter)")
    setup.ptol      = 1e-6
    setup.vtol      = 1e-7
    setup.btol      = 1e-6
    setup.nowarn    = 0
    setup.adj       = .
    b               = .z
}

void `MAIN'::clear()
{
    mu   = J(1,0,.)
    adj  = noadj = J(0,0,.)
    tau  = wsum = .
    b    = .z
    a    = .
    xb   = wbal = J(0,1,.)
    madj = J(1,0,.)
    loss = balanced = value = iter = conv = .
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
    setup.k = cols(X)
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
    // if scale is set by user
    if (setup.scale=="user") {
        if (setup.k!=length(scale)) _error(3200, "X not conformable with scale")
    }
    else scale = J(1,0,.) // (clear scale)
    // target moments
    if (setup.k!=cols(X0)) _error(3200, "X and Xref not conformable")
    setup.m0 = mean(X0, w0)
    // identify collinear terms
    setup.m  = mean(X, w)
    if (setup.k==0) {
        setup.omit = J(0,1,.)
        setup.s = J(1,0,.)
        setup.k_omit = 0
    }
    else {
        CP = quadcrossdev(X, setup.m, w, X, setup.m)
        setup.s = sqrt(diagonal(CP)' / setup.W)
        setup.omit = (diagonal(invsym(CP)):==0) // or: diagonal(invsym(CP, 1..setup.k)):==0
        setup.k_omit = sum(setup.omit)
    }
    // clear results
    setup.s0 = J(1,0,.) // (scale of refdata will be set later only if needed)
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

`RR' `MAIN'::s() return(setup.s)

`RR' `MAIN'::sref() {
    if (length(setup.s0)) return(setup.s0)
    if (nodata()) return(setup.s0)
    setup.s0 = sqrt(diagonal(quadcrossdev(*setup.X0, setup.m0, *setup.w0, 
        *setup.X0, setup.m0))'/setup.W0)
    return(setup.s0)
}

`Int' `MAIN'::k() return(setup.k)

`RC'  `MAIN'::omit() return(setup.omit)

`Int' `MAIN'::k_omit() return(setup.k_omit)

// settings -------------------------------------------------------------------

`T' `MAIN'::tau(| `TS' tau)
{
    if (args()==0) {
        if (tau<.) return(tau)
        if (nodata()) return(setup.tau)
        if      (setup.tau=="Wref") tau = setup.W0
        else if (setup.tau=="W")    tau = setup.W
        else if (setup.tau=="Nref") tau = setup.N0
        else if (setup.tau=="N")    tau = setup.N
        else                        tau = setup.tau
        return(tau)
    }
    if (setup.tau==tau) return // no change
    if (isstring(tau)) {
        if (!anyof(("Wref", "W", "Nref", "N"), tau)) {
            printf("{err}'%s' not allowed\n", tau)
            _error(3498)
        }
    }
    else if (tau<=0 | tau>=.) _error(3498, "setting out of range")
    setup.tau = tau
    clear()
}

`T' `MAIN'::scale(| `T' scale0)
{
    `RR' s
    
    if (args()==0) {
        if (length(scale)) return(scale)
        if (setup.scale=="user") return(scale)
        if (nodata()) return(setup.scale)
        if      (setup.scale=="main") scale = s()
        else if (setup.scale=="ref")  scale = sref()
        else if (setup.scale=="avg")  scale = (s() + sref()) / 2
        else if (setup.scale=="wavg") scale = (s()*W() + sref()*Wref()) / 
                                              (W() + Wref())
        else if (setup.scale=="pooled") {
            scale = sqrt(diagonal(mm_variance0(X() \ Xref(), 
                (rows(w())==1  ? J(N(), 1, w())   : w()) \ 
                (rows(wref())==1 ? J(Nref(), 1, wref()) : wref())))')
        }
        _editvalue(scale, 0, 1)
        return(scale)
    }
    if (isstring(scale0)) {
        if (length(scale0)!=1) _error(3200)
        if (setup.scale==scale0) return // no change
        if (!anyof(("main", "ref", "avg", "wavg","pooled"), scale0)) {
            printf("{err}'%s' not allowed\n", scale0)
            _error(3498)
        }
        setup.scale = scale0
        scale = J(1,0,.)
    }
    else {
        s = editvalue(vec(scale0)', 0, 1)
        if (scale==s) return // no change
        if (missing(s)) _error(3351)
        if (any(s:<=0)) _error(3498, "scale must be positive")
        if (nodata()==0) {
            if (setup.k!=length(s)) _error(3200, "scale not conformable with X")
        }
        setup.scale = "user"
        scale = s
    }
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

`T' `MAIN'::ltype(| `SS' ltype)
{
    if (args()==0) return(setup.ltype)
    if (setup.ltype==ltype) return // no change
    if (!anyof(("reldif", "absdif", "norm"), ltype)) {
        printf("{err}'%s' not allowed\n", ltype)
        _error(3498)
    }
    setup.ltype = ltype
    clear()
}

`T' `MAIN'::etype(| `SS' etype)
{
    if (args()==0) return(setup.etype)
    if (setup.etype==etype) return // no change
    if (!anyof(("bl","wl","mm","mma"), etype)) {
        printf("{err}'%g' not allowed\n", etype)
        _error(3498)
    }
    setup.etype = etype
    clear()
}

`T' `MAIN'::alteval(| `Bool' alteval) // for backward compatibility
{
    if (args()==0) return(setup.etype=="wl")
    if (setup.etype==(alteval!=0 ? "wl" : "bl")) return // no change
    setup.etype = (alteval!=0 ? "wl" : "bl")
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

`T' `MAIN'::difficult(| `Bool' difficult)
{
    if (args()==0) return(setup.difficult)
    if (setup.difficult==(difficult!=0)) return // no change
    setup.difficult = (difficult!=0)
    clear()
}

`T' `MAIN'::nostd(| `Bool' nostd)
{
    if (args()==0) return(setup.nostd)
    if (setup.nostd==(nostd!=0)) return // no change
    setup.nostd = (nostd!=0)
    clear()
}

`T' `MAIN'::nowarn(| `Bool' nowarn)
{
    if (args()==0) return(setup.nowarn)
    if (setup.nowarn==(nowarn!=0)) return // no change
    setup.nowarn = (nowarn!=0)
    clear()
}

// target means ---------------------------------------------------------------

`T' `MAIN'::adj(| `RM' adj0)
{
    `RR' adj
    
    if (args()==0) {
        if (nodata()) return(setup.adj)
        _setadj()
        return(this.adj)
    }
    adj = mm_unique(trunc(vec(adj0)))'
    if (length(setup.noadj)==0) {
        if (setup.adj==adj) return // no change
    }
    if (adj!=.) {
        if (any(adj:<1)) _error(3498, "setting out of range")
    }
    setup.adj   = adj
    setup.noadj = J(1,0,.)
    clear()
}

`T' `MAIN'::noadj(| `RM' noadj0)
{
    `RR' noadj
    
    if (args()==0) {
        if (nodata()) return(setup.noadj)
        _setadj()
        return(this.noadj)
    }
    noadj = mm_unique(trunc(vec(noadj0)))'
    if (setup.adj==.) {
        if (setup.noadj==noadj) return // no change
    }
    if (any(noadj:<1)) _error(3498, "setting out of range")
    setup.adj   = .
    setup.noadj = noadj
    clear()
}

void `MAIN'::_setadj() // assumes that data has been set
{
    `Int'  i, j, k
    `IntR' p
    
    if (rows(adj)) return // already set
    // zero variables
    k = k()
    if (k==0) {
        adj = noadj = J(1,0,.)
        return
    }
    // case 1: noadj() has been set
    if (length(setup.noadj)) {
        p = J(1,k,1)
        for (i=length(setup.noadj); i; i--) {
            j = setup.noadj[i]
            if (j>k) continue // be tolerant and ignore invalid subscripts
            p[j] = 0
        }
    }
    // case 2: adj() has been set
    else if (setup.adj!=.) {
        p = J(1,k,0)
        for (i=length(setup.adj); i; i--) {
            j = setup.adj[i]
            if (j>k) continue // be tolerant and ignore invalid subscripts
            p[j] = 1
        }
    }
    // case 3: default (adjust all)
    else p = 1
    // fill in adj/noadj
    if (allof(p, 1)) {
        adj   = .
        noadj = J(1,0,.)
    }
    else {
        adj   = select(1..k, p)
        noadj = select(1..k,!p)
    }
}

`RR' `MAIN'::mu()
{
    if (length(mu)) return(mu)
    if (nodata())   return(mu)
    if (adj()==.) {
        mu = setup.m0
    }
    else {
        mu = setup.m
        if (length(adj())) mu[adj()] = setup.m0[adj()]
    }
    return(mu)
}

// results --------------------------------------------------------------------

`RC' `MAIN'::b()
{
    if (b==.z) Fit()
    return(b)
}

`RS' `MAIN'::a()
{
    if (b==.z) Fit()
    return(a)
}

`RC' `MAIN'::wbal()
{
    if (b==.z) Fit()
    return(wbal)
}

`RC' `MAIN'::xb()
{
    if (b==.z) Fit()
    return(xb :+ a)
}

`RC' `MAIN'::pr()
{
    if (b==.z) Fit()
    return(invlogit(xb :+ (a + ln(Wref()/tau()))))
}

`RR' `MAIN'::madj()
{
    if (length(madj)) return(madj)
    if (b==.z) Fit()
    madj = mean(X(), wbal)
    return(madj)
}

`RS' `MAIN'::wsum()
{
    if (wsum<.) return(wsum)
    if (b==.z) Fit()
    wsum = quadsum(wbal)
    return(wsum)
}

`RS' `MAIN'::loss()
{
    if (b==.z) Fit()
    return(loss)
}

`RS' `MAIN'::balanced()
{
    if (b==.z) Fit()
    return(balanced)
}

`RS' `MAIN'::value()
{
    if (b==.z) Fit()
    return(value)
}

`Int' `MAIN'::iter()
{
    if (b==.z) Fit()
    return(iter)
}

`Bool' `MAIN'::converged()
{
    if (b==.z) Fit()
    return(conv)
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
    `Bool' nofit
    
    // optimize
    if (nodata()) _error(3498, "data not set")
    nofit = 0
    if ((k()-k_omit())<=0)   nofit = 1     // no covariates
    else if (!length(adj())) nofit = 1     // no adjustments
    else if (length(noadj()) & k_omit()) { // no adjustments among non-omitted
        nofit = all(omit()[adj()])
    }
    if (nofit) {
        b = J(k(), 1, 0)
        iter = value = 0
        conv = 1
    }
    else _Fit_b() // fit coefficients
    _Fit_a()      // compute intercept (and balancing weights)
    
    // check balancing
    if (etype()!="bl" | nostd()==0 | k_omit()) {
        // compute balancing loss using raw data
        loss = _mm_ebalance_loss(ltype(), mean(X():-mu(), wbal), mu())
    }
    else loss = value
    if (trace()!="none") {
        printf("{txt}Final fit:     balancing loss = {res}%10.0g\n", loss)
    }
    balanced = (loss<btol())
    if (balanced) return
    if (nowarn()) return
    display("{err}balance not achieved")
}

void `MAIN'::_Fit_a()
{
    `RS' ul
    
    xb = X() * b
    ul = max(xb) // set exp(max)=1 to avoid numerical overflow
    a  = ln(tau()) - ln(quadsum(w() :* exp(xb :- ul))) - ul
    wbal = w() :* exp(xb :+ a)
}

void `MAIN'::_Fit_b()
{
    `RR'   beta
    `IntC' p
    `T'    S
    
    // setup
    if (k_omit()) p = select(1::k(), omit():==0)
    S = optimize_init()
    optimize_init_which(S, "min")
    optimize_init_technique(S, "nr")
    optimize_init_evaluatortype(S, "d2")
    optimize_init_conv_ignorenrtol(S, "on")
    optimize_init_tracelevel(S, trace())
    optimize_init_conv_maxiter(S, maxiter())
    optimize_init_conv_ptol(S, ptol())
    optimize_init_conv_vtol(S, vtol())
    optimize_init_singularHmethod(S, difficult() ? "hybrid" : "")
    optimize_init_conv_warning(S, nowarn() ? "off" : "on")
    if (etype()=="mma") {
        optimize_init_evaluator(S, &_mm_ebalance_mma())   // gmm w/ alpha
        optimize_init_valueid(S, "criterion Q(p)")
        optimize_init_params(S, J(1, k()-k_omit()+1, 0)) // starting values
        optimize_init_argument(S, 1, _Fit_b_X(p))        // data
        optimize_init_argument(S, 2, _Fit_b_mu(p))       // target moments
        optimize_init_argument(S, 3, w()/W())            // norm. base weights
        optimize_init_argument(S, 4, tau()/W())          // normalized tau
    }
    else if (etype()=="mm") {
        optimize_init_evaluator(S, &_mm_ebalance_mm())   // gmm
        optimize_init_valueid(S, "criterion Q(p)")
        optimize_init_params(S, J(1, k()-k_omit(), 0))   // starting values
        optimize_init_argument(S, 1, _Fit_b_X(p))        // data
        optimize_init_argument(S, 2, _Fit_b_mu(p))       // target moments
        optimize_init_argument(S, 3, w())                // base weights
    }
    else if (etype()=="wl") {
        optimize_init_evaluator(S, &_mm_ebalance_lw())   // sum(ln(w))
        optimize_init_valueid(S, "criterion L(w)")
        optimize_init_params(S, J(1, k()-k_omit(), 0))   // starting values
        optimize_init_argument(S, 1, _Fit_b_Xc(p))       // centered data
        optimize_init_argument(S, 2, w())                // base weights
    }
    else {
        optimize_init_evaluator(S, &_mm_ebalance_bl())  // balance loss
        optimize_init_valueid(S, "balancing loss")
        optimize_init_params(S, J(1, k()-k_omit(), 0))  // starting values
        optimize_init_argument(S, 1, _Fit_b_Xc(p))      // centered data
        optimize_init_argument(S, 2, _Fit_b_mu(p))      // target moments
        optimize_init_argument(S, 3, w())               // base weights
        optimize_init_argument(S, 4, ltype())           // loss type
    }
    
    // run optimizer
    (void) _optimize(S)
    if (optimize_result_errorcode(S)) {
        errprintf("{p}\n")
        errprintf("%s\n", optimize_result_errortext(S))
        errprintf("{p_end}\n")
        exit(optimize_result_returncode(S))
    }
    
    // obtain results
    beta = optimize_result_params(S)
    if (etype()=="mma") beta = beta[|1\length(beta)-1|] // discard alpha
    if (k_omit()) {
        b = J(k(), 1, 0)
        if (nostd()) b[p] = beta'
        else         b[p] = (beta :/ scale()[p])'
    }
    else {
        if (nostd()) b = beta'
        else         b = (beta :/ scale())'
    }
    iter  = optimize_result_iterations(S)
    value = optimize_result_value(S)
    conv  = optimize_result_converged(S)
}

`RM' `MAIN'::_Fit_b_X(`IntC' p)
{
    if (k_omit()) {
        if (nostd()) return(X()[,p])
        return(X()[,p] :/ scale()[p])
    }
    if (nostd()) return(X())
    return(X() :/ scale())
}

`RM' `MAIN'::_Fit_b_Xc(`IntC' p)
{
    if (k_omit()) {
        if (nostd()) return(X()[,p] :- mu()[p])
        return((X()[,p] :- mu()[p]) :/ scale()[p])
    }
    if (nostd()) return(X() :- mu())
    return((X() :- mu()) :/ scale())
}

`RR' `MAIN'::_Fit_b_mu(`IntC' p)
{
    if (k_omit()) {
        if (nostd()) return(mu()[p])
        return((mu() :/ scale())[p])
    }
    if (nostd()) return(mu())
    return(mu() :/ scale())
}

`RS' _mm_ebalance_loss(`SS' ltype, `RR' d, `RR' mu)
{
    if (ltype=="absdif") return(max(abs(d)))
    if (ltype=="norm") return(sqrt(d*d'))
    return(mreldif(d+mu, mu)) // ltype=="reldif"
}

void _mm_ebalance_bl(`Int' todo, `RR' b, `RM' X, `RR' mu, `RC' w0, `SS' ltype,
    `RS' v, `RR' g, `RM' H)
{   // evaluator based on balance loss
    `RS' W
    `RC' w
    
    w = X * b'
    w = w0 :* exp(w :- max(w)) // avoid numerical overflow
    W = quadsum(w)
    g = quadcross(w, X) / W
    v = _mm_ebalance_loss(ltype, g, mu)
    if (todo==2) H = quadcross(X, w, X) / W
}

void _mm_ebalance_lw(`Int' todo, `RR' b, `RM' X, `RC' w0,
    `RS' v, `RR' g, `RM' H)
{   // evaluator using sum(ln(weights)) as criterion
    `RS' W
    `RC' w
    
    w = X * b'
    W = max(w)
    w = w0 :* exp(w :- W) // avoid numerical overflow
    v = ln(quadsum(w)) + W
    if (todo>=1) {
        W = quadsum(w)
        g = quadcross(w, X) / W
        if (todo==2) H = quadcross(X, w, X) / W
    }
}

void _mm_ebalance_mm(`Int' todo, `RR' b, `RM' X, `RR' mu, `RC' w0,
    `RS' v, `RR' g, `RM' H)
{   // gmm type evaluator
    `RC' w
    `RR' d
    `RM' h, G
    
    w = X * b'
    w = w0 :* exp(w :- max(w))   // avoid numerical overflow
    w = w :/ quadcolsum(w)
    h = w :* (X :- mu)
    d = quadcolsum(h)
    v = d * d'
    if (todo>=1) {
        G = quadcross(h, X :- quadcolsum(w:*X))
        g = d * G
        if (todo==2) H = G'G 
    }
}

void _mm_ebalance_mma(`Int' todo, `RR' b, `RM' X, `RR' mu, `RC' w0, `RS' tau, 
    `RS' v, `RR' g, `RM' H)
{   // gmm type evaluator including alpha in optimization problem
    `RS' d_a
    `RC' w, h_a
    `RR' d_b
    `RM' h_b, G
    
    w   = w0 :* exp(X*b[|1\cols(b)-1|]' :+ b[cols(b)])
    h_b = w :* (X :- mu)
    h_a = w :- w0:*tau  // = w0 * (e(xb+a) - tau/sum(w0))
    d_b = quadcolsum(h_b)
    d_a = quadsum(h_a)
    v   = (d_b, d_a) * (d_b, d_a)'
    if (todo>=1) {
        G = (quadcross(h_b, X), d_b') \ (quadcross(w, X), colsum(w))
        g = (d_b, d_a) * G
        if (todo==2) H = G'G 
    }
}

// influence functions --------------------------------------------------------

void `MAIN'::_IF_b()
{
    if (b==.z) Fit()
    _mm_ebalance_IF_b(IF, X(), Xref(), w(), wbal, madj(), mu(), tau(), W(),
        Wref(), adj(), noadj(), omit())
}

void `MAIN'::_IF_a()
{
    if (rows(IF.b)==0) _IF_b()
    _mm_ebalance_IF_a(IF, X(), w(), wbal, tau(), W())
}

void _mm_ebalance_IF_b(`If' IF, `RM' X, `RM' Xref, `RC' w, `RC' wbal, `RR' madj,
    `RR' mu, `RS' tau, `RS' W, `RS' Wref, `IntR' adj, `IntR' noadj,
    | `BoolC' omit)
{   // using "alternative approach" formulas
    `Int'  k
    `IntC' p
    `RM'   G, Q
    
    G = quadcrossdev(X, mu, wbal/tau, X, madj)
    if (length(omit)==0) {
        // no information on omitted terms; use qrinv()
        Q = -qrinv(G)
    }
    else if (any(omit)) {
        // discard omitted terms during inversion
        k = cols(X)
        p = select(1::k, omit:==0)
        Q = J(k, k, 0)
        if (length(p)) Q[p,p] = -luinv(G[p,p])
    }
    else {
        // no omitted terms; save to use luinv()
        Q = -luinv(G)
    }
    IF.b  = (wbal/tau):/w :* (X :- madj) * Q' // [sic!] using madj instead of mu
    if (adj==.) IF.b0 = (-1/Wref) * (Xref :- mu) * Q'
    else {
        IF.b0 = J(rows(Xref), cols(Xref), 0)
        if (length(adj)) IF.b0[, adj] = 
            (-1/Wref) * (Xref[,adj] :- mu[adj]) * Q[adj,adj]'
        if (length(noadj)) IF.b[,noadj] = 
            IF.b[,noadj] - (X[,noadj] :- mu[noadj])/W * Q[noadj,noadj]'
    }
}

void _mm_ebalance_IF_a(`If' IF, `RM' X, `RC' w, `RC' wbal, `RS' tau, `RS' W)
{    // using simplified IF assuming tau as fixed
    `RM' Q
    
    Q = quadcross(wbal, X)
    IF.a  = ((wbal:/w :- tau/W) + IF.b * Q') / -tau
    IF.a0 = (IF.b0 * Q') / -tau
}

end

