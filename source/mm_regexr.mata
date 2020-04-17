*! version 1.0.0  31jan2017  Ben Jann
version 9.2
mata:

string matrix mm_regexr(string matrix x, string matrix y, string matrix z)
{
    string matrix       res
    real scalar         r, R, c, C
    transmorphic scalar rx, cx, ry, cy, rz, cz
    pragma unset        r
    pragma unset        c

    R = max((rows(x),rows(y),rows(z)))
    C = max((cols(x),cols(y),cols(z)))
    rx = (rows(x)==1 ? &1 : (rows(x)<R ? _error(3200) : &r))
    cx = (cols(x)==1 ? &1 : (cols(x)<C ? _error(3200) : &c))
    ry = (rows(y)==1 ? &1 : (rows(y)<R ? _error(3200) : &r))
    cy = (cols(y)==1 ? &1 : (cols(y)<C ? _error(3200) : &c))
    rz = (rows(z)==1 ? &1 : (rows(z)<R ? _error(3200) : &r))
    cz = (cols(z)==1 ? &1 : (cols(z)<C ? _error(3200) : &c))
    res = J(R,C, "")
    for (r=1;r<=R;r++) {
        for (c=1;c<=C;c++) {
            res[r,c] = _mm_regexr(x[*rx,*cx], y[*ry,*cy], z[*rz,*cz])
        }
    }
    return(res)
}

string scalar _mm_regexr(string scalar s0, string scalar from, string scalar to)
{
    real scalar         i, j
    string scalar       s, t, BSLASH
    string rowvector    sub
    pragma unset        s
    
    if (regexm(s0, from)) {
        sub = regexs()
        sub = (sub, J(1, 10-cols(sub), ""))
        BSLASH = "\"
        for (i=1; i<=strlen(to); i++) {
            if ((t=substr(to,i,1))==BSLASH) {
                i++
                if ((t=substr(to,i,1))==BSLASH) s = s + BSLASH      // "\\"
                else if ((j=strtoreal(t))<.)    s = s + sub[j+1]    // "\#"
                else                            s = s + BSLASH + t
            }
            else s = s + t
        }
        return(subinstr(s0, sub[1], s, 1))
    }
    else return(s0)
}

end
