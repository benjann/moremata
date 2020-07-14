*! version 1.0.2  09jul2020  Ben Jann

version 10.0
mata:

numeric matrix _mm_colrunsum10(numeric matrix A, real scalar mis)
{
    //version 10
    real scalar i
    numeric matrix B

    if (cols(A)==1) return(runningsum(A, mis))

    if (isfleeting(A)) {
        for (i=1; i<=cols(A); i++) A[,i] = runningsum(A[,i], mis)
        return(A)
    }
    B = A
    for (i=1; i<=cols(B); i++) B[,i] = runningsum(B[,i], mis)
    return(B)
}

numeric matrix _mm_quadcolrunsum10(numeric matrix A, real scalar mis)
{
    //version 10
    real scalar i
    numeric matrix B

    if (cols(A)==1) return(quadrunningsum(A, mis))

    if (isfleeting(A)) {
        for (i=1; i<=cols(A); i++) A[,i] = quadrunningsum(A[,i], mis)
        return(A)
    }
    B = A
    for (i=1; i<=cols(B); i++) B[,i] = quadrunningsum(B[,i], mis)
    return(B)
}

end
