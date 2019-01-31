*! version 1.0.8  11jan2008  Ben Jann
version 9.0       
mata:             
                  
numeric matrix mm_colrunsum(numeric matrix A)
{
    numeric matrix B

    if (stataversion()>=1000) return(_mm_colrunsum_goto10(A))

    if (isfleeting(A)) {
        _mm_colrunsum(A)
        return(A)
    }
    _mm_colrunsum(B=A)
    return(B)
}

void _mm_colrunsum(numeric matrix Z)
{
    real scalar i

    _editmissing(Z, 0)
    for (i=2; i<=rows(Z); i++) Z[i,] = Z[i-1,] + Z[i,]

//    // running sum using mean-update formula; see Gould 2006, SJ 6(4)
//    _editmissing(Z, 0)
//    if (rows(Z)<2) return
//    for (i=2; i<=rows(Z); i++) Z[i,] = Z[i-1,] + (Z[i,]-Z[i-1,])/i
//    Z = Z :* (1::rows(Z))

}

numeric matrix _mm_colrunsum_goto10(numeric matrix Z)
{
    return(_mm_colrunsum10(Z))
}


end
