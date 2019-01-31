*! version 1.0.1, Ben Jann, 14apr2006
version 9.1
mata:

real matrix mm_rbinomial(real matrix n, real matrix p)
{
	real matrix x
	real scalar r, R, c, C
	transmorphic scalar rn, cn, rp, cp

	R = max((rows(n),rows(p)))
	C = max((cols(n),cols(p)))
	rn = (rows(n)==1 ? &1 : (rows(n)<R ? _error(3200) : &r))
	cn = (cols(n)==1 ? &1 : (cols(n)<C ? _error(3200) : &c))
	rp = (rows(p)==1 ? &1 : (rows(p)<R ? _error(3200) : &r))
	cp = (cols(p)==1 ? &1 : (cols(p)<C ? _error(3200) : &c))
	x = J(R,C,.)
	for (r=1;r<=R;r++) {
		for (c=1;c<=C;c++) {
			x[r,c] = _mm_rbinomial(n[*rn,*cn], p[*rp,*cp])
		}
	}
	return(x)
}

real scalar _mm_rbinomial(real scalar n, real scalar p)
{
	if (p<0|p>1) return(.)
	if (n<=0|(n-trunc(n))!=0) return(.)
//rejection technique
	if (n<50|p>.03) {
		return(colsum(uniform(n,1):<p))
	}
//geometric distribution technique
	real scalar sum, x
	sum = 0
	for (x=1;x<=n;x++) {
		sum = sum + ceil((ln(uniform(1,1))/ln(1-p))-1)
		if (sum > n-x) break
	}
	return(x-1)
}

end
