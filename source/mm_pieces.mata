*! version 1.0.2, Ben Jann, 01jun2015
version 9.2
mata:

string rowvector mm_pieces(
    string scalar s,
    real scalar w,
    | real scalar nobreak,
      real scalar ascii)
{
    real scalar         i, j, k, n, l, b, nobr
    string scalar c
    string rowvector    res

    nobr = ( args()>2 ? nobreak : 0 )
    if (stataversion()>=1400 & !(args()>3 ? ascii : 0)) {
        return(_mm_pieces_goto14(s, w, nobr))
    }
    l = strlen(s)
    if (l<2 | w>=l) return(strtrim(s))
    res = J(1, _mm_npieces(s, w, nobr, (args()>3 ? ascii : 0)), "")
    j = k = n = 0
    b = 1
    for (i=1; i<=l; i++) {
        c = substr(s, i, 1)
        if (j<1) { // skip to first nonblank character
            if (c==" ") {
                b++
                continue
            }
        }
        j++
        if (i==l) res[++n] = strtrim(substr(s, b, .))
        else {
            if (c==" ") k = i
            if (j>=w) {
                if (k<1) {
                    if (nobr) continue
                    k = i
                }
                else {
                    if (substr(s, i+1, 1)==" ") k = i
                }
                res[++n] = strtrim(substr(s, b, k-b+1))
                j = i - k
                b = k + 1
                k = 0
            }
        }
    }
    return(res)
}

real matrix mm_npieces(
    string matrix S,
    real scalar w,
    | real scalar nobreak,
      real scalar ascii)
{
    real scalar  i, j
    real matrix  res

    res = J(rows(S),cols(S),1)
    if (w>=.) return(res)
    for (i=1; i<=rows(S); i++) {
        for (j=1; j<=cols(S); j++) {
            res[i,j] = _mm_npieces(S[i,j], w, (args()>2 ? nobreak : 0), 
                                              (args()>3 ? ascii : 0))
        }
    }
    return(res)
}

real scalar _mm_npieces(
    string scalar s,
    real scalar w,
    | real scalar nobreak,
      real scalar ascii)
{
    real scalar   i, j, k, n, l, nobr
    string scalar c

    nobr = ( args()>2 ? nobreak : 0 )
    if (stataversion()>=1400 & !(args()>3 ? ascii : 0)) {
        return(_mm_npieces_goto14(s, w, nobr))
    }
    l = strlen(s)
    if (l<2 | w>=l) return(1)
    j = k = n = 0
    for (i=1; i<=l; i++) {
        c = substr(s, i, 1)
        if (j<1) { // skip to first nonblank character
            if (c==" ") continue
        }
        j++
        if (i==l) n++
        else {
            if (c==" ") k = i
            if (j>=w) {
                if (k<1) {
                    if (nobr) continue
                    k = i
                }
                else {
                    if (substr(s, i+1, 1)==" ") k = i
                }
                n++
                j = i - k
                k = 0
            }
        }
    }
    if (n==0) n++
    return(n)
}

string rowvector _mm_pieces_goto14(string scalar s,
    real scalar w,
    real scalar nobr)
{
    return(_mm_pieces14(s, w, nobr))
}

real scalar _mm_npieces_goto14(string scalar s,
    real scalar w,
    real scalar nobr)
{
    return(_mm_npieces14(s, w, nobr))
}

end
