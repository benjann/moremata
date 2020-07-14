*! version 1.0.0  11jul2020  Ben Jann
version 9.2
mata:

real matrix mm_clip(real matrix X, real scalar min, real scalar max, | real scalar miss)
{
    real matrix R
    
    if (args()<4) miss = 0
    if (isfleeting(X)) {
        _mm_clip(X, min, max, miss)
        return(X)
    }
    R = X
    _mm_clip(R, min, max, miss)
    return(R)
}

void _mm_clip(real matrix X, real scalar min, real scalar max, | real scalar miss)
{
    real scalar i, j, r, x
    
    if (args()<4) miss = 0
    if (min>max) _error(3300)
    r = rows(X)
    if (miss) {
        for (j=cols(X); j; j--) {
            for (i=r; i; i--) {
                x = X[i,j]
                if (x>=.) continue
                if (x<min)      X[i,j] = min
                else if (x>max) X[i,j] = max
            }
        }
        return
    }
    for (j=cols(X); j; j--) {
        for (i=r; i; i--) {
            x = X[i,j]
            if (x<min)      X[i,j] = min
            else if (x>max) X[i,j] = max
        }
    }
}

real matrix mm_clipmin(real matrix X, real scalar min)
{
    real matrix R
    
    if (isfleeting(X)) {
        _mm_clipmin(X, min)
        return(X)
    }
    R = X
    _mm_clipmin(R, min)
    return(R)
}

void _mm_clipmin(real matrix X, real scalar min)
{
    real scalar i, j, r
    
    r = rows(X)
    for (j=cols(X); j; j--) {
        for (i=r; i; i--) {
            if (X[i,j]<min) X[i,j] = min
        }
    }
}

real matrix mm_clipmax(real matrix X, real scalar max, | real scalar miss)
{
    real matrix R
    
    if (args()<3) miss = 0
    if (isfleeting(X)) {
        _mm_clipmax(X, max, miss)
        return(X)
    }
    R = X
    _mm_clipmax(R, max, miss)
    return(R)
}

void _mm_clipmax(real matrix X, real scalar max, | real scalar miss)
{
    real scalar i, j, r, x
    
    if (args()<3) miss = 0
    r = rows(X)
    if (miss) {
        for (j=cols(X); j; j--) {
            for (i=r; i; i--) {
                x = X[i,j]
                if (x>=.) continue
                if (x>max) X[i,j] = max
            }
        }
        return
    }
    for (j=cols(X); j; j--) {
        for (i=r; i; i--) {
            if (X[i,j]>max) X[i,j] = max
        }
    }
}

end
