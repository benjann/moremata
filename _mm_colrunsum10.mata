*! version 1.0.1  11jan2008  Ben Jann

version 10.0
mata:

numeric matrix _mm_colrunsum10(numeric matrix A)
{
    //version 10
    real scalar i
    numeric matrix B

    if (cols(A)==1) return(runningsum(A))

    if (isfleeting(A)) {
        for (i=1; i<=cols(A); i++) A[,i] = runningsum(A[,i])
        return(A)
    }
    B = A
    for (i=1; i<=cols(B); i++) B[,i] = runningsum(B[,i])
    return(B)
}

end
