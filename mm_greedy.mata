*! version 1.0.0  02may2019  Ben Jann
version 9.2
local Int   real scalar
local IntC  real colvector
local IntM  real matrix
local RS    real scalar
local RC    real colvector
local T     transmorphic
local TM    transmorphic matrix
local PidC  pointer(`Int') colvector
local PidM  pointer(`Int') matrix
local dist  pointer(function) scalar
mata:

`IntM' mm_greedy(`TM' T, `TM' C, `Int' n, `RS' calip, `dist' f, | `T' fopts)
{
    `RC'   d
    `PidM' ij
    pragma unset d
    
    // check input
    if (cols(T)!=cols(C)) _error(3200, "T and C must have same number of columns")
    if (cols(T)==0) _error(3200, "T and C must have at least one column")
    // make ID lookup table and compute distances
    if (calip>=.) ij = _mm_greedy_dist(d, T, C, f, fopts)
    else          ij = _mm_greedy_dist_calip(d, T, C, f, fopts, calip)
    // match and return IDs of matched controls
    if (n<=1 | n>=.) return(_mm_greedy_match_1(ij, d, rows(T)))
    return(_mm_greedy_match_n(ij, d, rows(T), n))
}
`PidM' _mm_greedy_dist(`RC' d, `TM' T, `TM' C, `dist' f, `T' fopts)
{
    `Int'  i, a, b, nT, nC, N
    `PidC' Ti, Ci
    `PidM' ij
    
    nT = rows(T); nC = rows(C)   // number of treated and controls
    N = nT * nC                  // number of combinations
    Ti = _mm_greedy_pid(nT)      // make ID pointers for treated
    Ci = _mm_greedy_pid(nC)      // make ID pointers for controls
    d = J(N, 1, .)               // prepare vector of distances
    ij = J(N, 2, NULL)           // prepare ID lookup table
    if (N==0) return(ij)         // T or C is void
    b = 0
    for (i=nT; i; i--) {         // loop over treated => compare to all controls
        a = b + 1                // start index of current batch
        b = b + nC               // end index of current batch
        d[|a\b|] = (*f)(T[i,], C, fopts)   // add distances to distance vector
        ij[|a,1\b,.|] = J(nC,1,Ti[i]), Ci  // add ID pointers to lookup table
    }
    return(ij)
}
`PidM' _mm_greedy_dist_calip(`RC' d, `TM' T, `TM' C, `dist' f, `T' fopts,
    `RS' calip)
{
    `Int'  i, a, b, nT, nC, N, n
    `PidC' Ti, Ci
    `PidM' ij
    `RC'   D
    `IntC' P, p
    
    nT = rows(T); nC = rows(C)    // number of treated and controls
    N = nT * nC                   // number of combinations
    Ti = _mm_greedy_pid(nT)       // make ID pointers for treated
    Ci = _mm_greedy_pid(nC)       // make ID pointers for controls
    d = J(N, 1, .)                // prepare vector of distances
    ij = J(N, 2, NULL)            // prepare ID lookup table
    if (N==0) return(ij)          // T or C is void
    b = 0; P = (1::nC)
    for (i=nT; i; i--) {          // loop over treated => compare to all controls
        D = (*f)(T[i,], C, fopts) // compute distances
        p = select(P, D:<=calip)  // index of distances within caliper
        n = length(p)             // number of distances within caliper
        if (n==0) continue        // no valid distances -> skip case
        a = b + 1                 // start index of current batch
        b = b + n                 // end index of current batch
        if (n==nC) {              // all distances are valid
            d[|a\b|] = D          // add distances to distance vector
            ij[|a,1\b,.|] = J(nC,1,Ti[i]), Ci // add ID pointers to lookup table
            continue
        }
        d[|a\b|] = D[p]          // add distances to distance vector
        ij[|a,1\b,.|] = J(n,1,Ti[i]), Ci[p] // add ID pointers to lookup table
    }
    if (b==N) return(ij)         // all distances valid; all rows filled
    if (b==0) {                  // no valid distances
        d = J(0, 1, .)
        return(J(0, 2, NULL))
    }
    d = d[|1\b|]                 // remove unused rows from d
    return(ij[|1,1\b,.|])        // remove unused rows from ij
}

`PidC' _mm_greedy_pid(`Int' n)
{
    `Int'  i
    `PidC' P
    
    P = J(n, 1, NULL)                 // prepare container for ID pointers
    for (i=n; i; i--) P[i] = &(i*1)   // generate ID pointers
    return(P)
}

`IntC' _mm_greedy_match_1(`PidM' ij, `RC' d, `Int' nT)
{
    `Int'  i, j, ii, jj, N, L
    `IntC' p, M
    
    M = J(nT, 1, .)                        // prepare container for matched IDs
    L = nT                                 // maximum number of matches
    N = rows(d)                            // number of distances
    p = order(d, 1)                        // obtain sort order of distances
    for (i=1; i<=N; i++) {
        j = p[i]                           // get next position
        if ((ii = *(ij[j,1]))==.) continue // treatment case already matched
        if ((jj = *(ij[j,2]))==.) continue // control case already used
        M[ii] = jj                         // copy matched control ID
        *(ij[j,1]) = .                     // delete treatment ID
        *(ij[j,2]) = .                     // delete control ID
        if (!(--L)) break                  // maximum reached; no treatment IDs left
    }
    return(M)
}

`IntM' _mm_greedy_match_n(`PidM' ij, `RC' d, `Int' nT, `Int' n)
{
    `Int'  i, j, ii, jj, N, L, nc
    `IntC' p, NC
    `IntM' M
    
    NC = J(nT, 1, 0)                       // counter for number of matches
    M = J(nT, n, .)                        // container for matched IDs
    L = nT * n                             // maximum number of matches
    N = rows(d)                            // number of distances
    p = order(d, 1)                        // obtain sort order of distances
    for (i=1; i<=N; i++) {
        j = p[i]                           // get next position
        if ((ii = *(ij[j,1]))==.) continue // treatment case already matched
        if ((jj = *(ij[j,2]))==.) continue // control case already used
        nc = NC[ii]; nc++                  // update number of matches
        M[ii, nc] = jj                     // copy matched control ID
        if (nc==n) *(ij[j,1]) = .          // delete treatment ID
        *(ij[j,2]) = .                     // delete control ID
        NC[ii] = nc                        // store number of matches
        if (!(--L)) break                  // maximum reached; no treatment IDs left
    }
    return(M)
}

end


