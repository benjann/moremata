*! version 1.0.0  01may2025  Ben Jann

version 9.2
mata:

string matrix mm_read_csv(string scalar filename, | string scalar del,
    real scalar notrim, real scalar nostrip, real scalar l1, real scalar l2)
{
    real scalar      a, b, c, ci, i, j, l, q
    string scalar    s, dq
    string rowvector row
    string matrix    S
    
    if (args()<2) del = ","
    if (args()<3) notrim = 0
    if (args()<4) nostrip = 0
    if (args()<5) l1  = 1
    dq = `"""'
    S = cat(filename, l1, l2)
    i = rows(S)
    if (!i) return(J(0,0,"")) // input is empty
    l = strlen(del)
    if (!l) { // del is empty string; return one column
        if (!notrim)  S = strtrim(S)
        if (!nostrip) _mm_read_csv_strip(S, dq)
        return(S)
    }
    c = ci = 1
    for (;i;i--) {
        row = J(1,ci,"")
        s = S[i,1]
        b = _mm_read_csv_pos(s, del, 1, 0) // first delimiter
        q = _mm_read_csv_pos(s, dq, 1, 0)  // first quote
        a = 1; j = 0
        while (1) {
            if (q<b) {
                q = _mm_read_csv_pos(s, dq, 1, q) // matching quote
                if (q<.) {
                    b = _mm_read_csv_pos(s, del, 1, q) // next delimiter
                    q = _mm_read_csv_pos(s, dq, 1, q)  // next quote
                    continue
                }
            }
            if (b<.) {
                _mm_read_csv_put(row, ci, j, s, a, b) // add token to row
                a = b + l
                b = _mm_read_csv_pos(s, del, l, b) // next delimiter
                continue
            }
            _mm_read_csv_put(row, ci, j, s, a, b) // add row to result
            if (j>c) {
                S = S, J(rows(S),j-c,"")
                c = j
                S[i,] = row
            }
            else S[i,] = row
            break
        }
    }
    if (!notrim)  S = strtrim(S)
    if (!nostrip) _mm_read_csv_strip(S, dq)
    return(S)
}

real scalar _mm_read_csv_pos(string scalar s, string scalar t, real scalar l,
    real scalar o)
{
    real scalar p
    
    if (o) {
        p = strpos(substr(s,o+l,.), t)
        if (p) return(p+o)
        return(.)
    }
    p = strpos(s, t)
    if (p) return(p)
    return(.)
}

void _mm_read_csv_put(string rowvector row, real scalar c, real scalar j,
    string scalar s, real scalar a, real scalar b)
{
    if (++j>c) {
        row = row, substr(s,a,b-a)
        c++
    }
    else row[j] = substr(s,a,b-a)
}

void _mm_read_csv_strip(string matrix S, string scalar dq)
{
    real scalar   i, j, r
    string scalar s
    
    r = rows(S)
    for (j=cols(S);j;j--) {
        for (i=r;i;i--) {
            s = S[i,j]
            if (substr(S[i,j],1,1)==dq) {
                if (substr(s,-1,1)==dq) {
                    // strip outer double quotes
                    S[i,j] = substr(s,2,strlen(s)-2)
                }
            }
        }
    }
}

end
