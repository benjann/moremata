*! version 1.0.2  23apr2021  Ben Jann
version 9.2

local opts "o1, o2, o3, o4, o5, o6, o7, o8, o9, o10"

mata:

real matrix mm_collapse(
    real matrix    X,              // data
    real colvector w,              // weights (or 1)
    real colvector ID,             // subgroup IDs
    |  pointer(function) scalar f, // function; default: mean()
      `opts')                      // arguments to pass through to f()
{
    real colvector p
    pointer scalar Xs, IDs, ws

    if ((rows(X)!=rows(ID)) | (rows(w)!=1 & rows(X)!=rows(w))) _error(3200)
    p     = order(ID, 1)
    if (isfleeting(X))      Xs  = &(X = X[p,])
    else                    Xs  = &(X[p,])
    if (rows(w)==1)         ws  = &w
    else if (isfleeting(w)) ws  = &(w = w[p])
    else                    ws  = &(w[p])
    if (isfleeting(ID))     IDs = &(ID = ID[p])
    else                    IDs = &(ID[p])

    if (args()==3)  return(_mm_collapse(*Xs, *ws, *IDs))
    if (args()==4)  return(_mm_collapse(*Xs, *ws, *IDs, f))
    if (args()==5)  return(_mm_collapse(*Xs, *ws, *IDs, f, o1))
    if (args()==6)  return(_mm_collapse(*Xs, *ws, *IDs, f, o1, o2))
    if (args()==7)  return(_mm_collapse(*Xs, *ws, *IDs, f, o1, o2, o3))
    if (args()==8)  return(_mm_collapse(*Xs, *ws, *IDs, f, o1, o2, o3, o4))
    if (args()==9)  return(_mm_collapse(*Xs, *ws, *IDs, f, o1, o2, o3, o4, o5))
    if (args()==10) return(_mm_collapse(*Xs, *ws, *IDs, f, o1, o2, o3, o4, o5, o6))
    if (args()==11) return(_mm_collapse(*Xs, *ws, *IDs, f, o1, o2, o3, o4, o5, o6, o7))
    if (args()==12) return(_mm_collapse(*Xs, *ws, *IDs, f, o1, o2, o3, o4, o5, o6, o7, o8))
    if (args()==13) return(_mm_collapse(*Xs, *ws, *IDs, f, o1, o2, o3, o4, o5, o6, o7, o8, o9))
                    return(_mm_collapse(*Xs, *ws, *IDs, f, `opts'))
}


real matrix _mm_collapse(
    real matrix    X,               // data
    real colvector w,               // weights (or 1)
    real colvector ID,              // subgroup IDs
    |  pointer(function) scalar f0, // function; default: mean()
      `opts')                       // arguments to pass through to f()
{
    real scalar                 i, j, c, a, b, ww
    real matrix                 info, res, key
    pointer(function) scalar    f
    transmorphic                setup

    if (rows(X)!=rows(ID))     _error(3200)
    ww      = rows(w)!=1
    if (ww & rows(w)!=rows(X)) _error(3200)
    if (rows(X)<1)             return(J(0, cols(X)+1, missingof(X)))

    f     = args()<4 ? &mean() : f0
    setup = mm_callf_setup(f, args()-4, `opts')

    info = _mm_panels(ID)
    key  = J(rows(info), 1, missingof(ID))
    res  = J(rows(info), c=cols(X), missingof(X))
    b = 0
    for (i=1; i<=rows(info); i++) {
        a = b + 1
        b = b + info[i]
        key[i] = ID[a]
        for (j=1; j<=c; j++) {
            res[i, j] =
                mm_callf(setup, X[|a,j \ b,j|], ww ? w[|a \ b|] : w)
        }
    }
    return(key, res)
}

real matrix mm_collapse2(
    real matrix    X,              // data
    real colvector w,              // weights (or 1)
    real colvector ID,             // subgroup IDs
    |  pointer(function) scalar f, // function; default: mean()
      `opts')                      // arguments to pass through to f()
{
    real colvector p
    pointer scalar Xs, IDs, ws
    real matrix    res

    if ((rows(X)!=rows(ID)) | (rows(w)!=1 & rows(X)!=rows(w))) _error(3200)
    p     = order(ID, 1)
    if (isfleeting(X))      Xs  = &(X = X[p,])
    else                    Xs  = &(X[p,])
    if (rows(w)==1)         ws  = &w
    else if (isfleeting(w)) ws  = &(w = w[p])
    else                    ws  = &(w[p])
    if (isfleeting(ID))     IDs = &(ID = ID[p])
    else                    IDs = &(ID[p])

    res = J(rows(X), cols(X), missingof(X))
    if      (args()==3)  res[p,] = _mm_collapse2(*Xs, *ws, *IDs)
    else if (args()==4)  res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f)
    else if (args()==5)  res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f, o1)
    else if (args()==6)  res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f, o1, o2)
    else if (args()==7)  res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f, o1, o2, o3)
    else if (args()==8)  res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f, o1, o2, o3, o4)
    else if (args()==9)  res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f, o1, o2, o3, o4, o5)
    else if (args()==10) res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f, o1, o2, o3, o4, o5, o6)
    else if (args()==11) res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f, o1, o2, o3, o4, o5, o6, o7)
    else if (args()==12) res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f, o1, o2, o3, o4, o5, o6, o7, o8)
    else if (args()==13) res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f, o1, o2, o3, o4, o5, o6, o7, o8, o9)
    else                 res[p,] = _mm_collapse2(*Xs, *ws, *IDs, f, `opts')
    return(res)
}

real matrix _mm_collapse2(
    real matrix    X,               // data
    real colvector w,               // weights (or 1)
    real colvector ID,              // subgroup IDs
    |  pointer(function) scalar f0, // function; default: mean()
      `opts')                       // arguments to pass through to f()
{
    real scalar                 i, j, c, a, b, ww
    real matrix                 info, res
    pointer(function) scalar    f
    transmorphic                setup

    if (rows(X)!=rows(ID))     _error(3200)
    ww      = rows(w)!=1
    if (ww & rows(w)!=rows(X)) _error(3200)
    if (rows(X)<1)             return(J(0, cols(X), missingof(X)))

    f     = args()<4 ? &mean() : f0
    setup = mm_callf_setup(f, args()-4, `opts')

    info = _mm_panels(ID)
    res  = J(rows(X), c=cols(X), missingof(X))
    b = 0
    for (i=1; i<=rows(info); i++) {
        a = b + 1
        b = b + info[i]
        for (j=1; j<=c; j++) {
            res[|a,j \ b,j|] = J(info[i], 1,
                mm_callf(setup, X[|a,j \ b,j|], ww ? w[|a \ b|] : w))
        }
    }
    return(res)
}

end
