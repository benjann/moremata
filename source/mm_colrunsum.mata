*! version 1.0.9  09jul2020  Ben Jann
version 9.0
mata:

numeric matrix mm_colrunsum(numeric matrix A, | real scalar mis, real scalar qd)
{
    numeric matrix B

    if (args()<2) mis = 0
    if (stataversion()>=1000) {
        if (args()<3) qd = 0
        return(_mm_colrunsum_goto10(A, mis, qd))
    }
    if (isfleeting(A)) {
        _mm_colrunsum(A, mis)
        return(A)
    }
    _mm_colrunsum(B=A, mis)
    return(B)
}

void _mm_colrunsum(numeric matrix Z, real scalar mis)
{
    real scalar i

    if (mis==0) _editmissing(Z, 0)
    for (i=2; i<=rows(Z); i++) Z[i,] = Z[i-1,] + Z[i,]

//    // running sum using mean-update formula; see Gould 2006, SJ 6(4)
//    if (mis==0) _editmissing(Z, 0)
//    if (rows(Z)<2) return
//    for (i=2; i<=rows(Z); i++) Z[i,] = Z[i-1,] + (Z[i,]-Z[i-1,])/i
//    Z = Z :* (1::rows(Z))

}

numeric matrix _mm_colrunsum_goto10(numeric matrix Z, real scalar mis, real scalar qd)
{
    if (qd) return(_mm_quadcolrunsum10(Z, mis))
    return(_mm_colrunsum10(Z, mis))
}


end
