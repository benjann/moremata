*! version 1.0.1, Ben Jann, 14apr2006
version 9.1
mata:

real matrix mm_cebinomial(real matrix n, real matrix k, real matrix p)
{
	real matrix x
	real scalar r, R, c, C
	transmorphic scalar rn, cn, rk, ck, rp, cp

	R = max((rows(n),rows(k),rows(p)))
	C = max((cols(n),cols(k),cols(p)))
	rn = (rows(n)==1 ? &1 : (rows(n)<R ? _error(3200) : &r))
	cn = (cols(n)==1 ? &1 : (cols(n)<C ? _error(3200) : &c))
	rk = (rows(k)==1 ? &1 : (rows(k)<R ? _error(3200) : &r))
	ck = (cols(k)==1 ? &1 : (cols(k)<C ? _error(3200) : &c))
	rp = (rows(p)==1 ? &1 : (rows(p)<R ? _error(3200) : &r))
	cp = (cols(p)==1 ? &1 : (cols(p)<C ? _error(3200) : &c))
	x = J(R,C,.)
	for (r=1;r<=R;r++) {
		for (c=1;c<=C;c++) {
			x[r,c] = _mm_cebinomial(n[*rn,*cn], k[*rk,*ck], p[*rp,*cp])
		}
	}
	return(x)
}

real scalar _mm_cebinomial(real scalar n, real scalar k, real scalar p)
{
	real scalar e, e0, i

	if (p<0|p>1) return(.)                 //success prob. out of range
	if (n<=0|(n-trunc(n))!=0) return(.)    //n out of range or non-integer
	if (k<0|k>n|(k-trunc(k))!=0) return(.) //k out of range or non-integer
	if (k==0) return(n*p)
	e = 0
	e0 = 0
	for (i=1;i<=n-k;i++) {
		e = e + Binomial(n,k+i,p)
		if (e==e0) break
		e0 = e
	}
	return(e / Binomial(n,k,p) + k)
}
end
