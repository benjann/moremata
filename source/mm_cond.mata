*! version 1.0.0  29feb2008  Ben Jann
version 9.2
mata:

transmorphic matrix mm_cond(real matrix x, transmorphic matrix y, transmorphic matrix z)
{
        transmorphic matrix res
        real scalar r, R, c, C
        transmorphic scalar rx, cx, ry, cy, rz, cz

        if (eltype(y) != eltype(z)) _error(3250)
        R = max((rows(x),rows(y),rows(z)))
        C = max((cols(x),cols(y),cols(z)))
        rx = (rows(x)==1 ? &1 : (rows(x)<R ? _error(3200) : &r))
        cx = (cols(x)==1 ? &1 : (cols(x)<C ? _error(3200) : &c))
        ry = (rows(y)==1 ? &1 : (rows(y)<R ? _error(3200) : &r))
        cy = (cols(y)==1 ? &1 : (cols(y)<C ? _error(3200) : &c))
        rz = (rows(z)==1 ? &1 : (rows(z)<R ? _error(3200) : &r))
        cz = (cols(z)==1 ? &1 : (cols(z)<C ? _error(3200) : &c))
        res = J(R,C, missingof(y))
        for (r=1;r<=R;r++) {
                for (c=1;c<=C;c++) {
                        res[r,c] = (x[*rx,*cx] ? y[*ry,*cy] : z[*rz,*cz])
                }
        }
        return(res)
}

end
