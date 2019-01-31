*! version 1.0.0, Ben Jann, 07jun2006
version 9.2
mata:

transmorphic matrix mm_expand(transmorphic matrix X, real vector wr,
 | real vector wc, real scalar sort)
{
	transmorphic matrix Xnew

	if (args()<3) wc = 1
	if (args()<4) sort = 0
	if (isfleeting(X)) {
		_mm_expand(X, wr, wc, sort)
		return(X)
	}
	else {
		_mm_expand(Xnew=X, wr, wc, sort)
		return(Xnew)
	}
}

void _mm_expand(transmorphic matrix X, real vector wr,
 real vector wc, real scalar sort)
{
	real scalar i, b, e, n, add
	real colvector f
	pointer scalar fi

// expand rows
	if (length(wr)==1 & sort==0) _mm_repeat(X, wr, 1)
	else {
		f = trunc(wr:-1) :* (wr:>0)
		_editmissing(f,0)
		n = rows(X)
		if (length(f)==1) {
			fi = &1
			add = n*f
		}
		else {
			if (n!=length(f)) _error(3200)
			fi = &i
			add = sum(f)
		}
		if (add>0) {
			X = X \ J(add,cols(X),missingof(X))
			if (cols(X)>0) {
				if (sort) f = f :+ 1
				b = rows(X)+1
				for (i=n; i>=1; i--) {
					if (f[*fi]<1) continue
					e = b - 1
					b = b - f[*fi]
					X[|b,1 \ e,.|] = X[J(f[*fi],1,i), .]
				}
			}
		}
	}
// expand columns
	if (length(wc)==1 & sort==0) {
		_mm_repeat(X, 1, wc)
		return
	}
	f = trunc(wc:-1) :* (wc:>0)
	_editmissing(f,0)
	n = cols(X)
	if (length(f)==1) {
		fi = &1
		add = n*f
	}
	else {
		if (n!=length(f)) _error(3200)
		fi = &i
		add = sum(f)
	}
	if (add>0) {
		X = X , J(rows(X),add,missingof(X))
		if (rows(X)>0) {
			if (sort) f = f :+ 1
			b = cols(X)+1
			for (i=n; i>=1; i--) {
				if (f[*fi]<1) continue
				e = b - 1
				b = b - f[*fi]
				X[|1,b \ .,e|] = X[., J(1,f[*fi],i)]
			}
		}
	}
}

transmorphic matrix mm_repeat(transmorphic matrix X, real scalar wr,
 | real scalar wc)
{
	transmorphic matrix Xnew

	if (args()<3) wc = 1
	if (isfleeting(X)) {
		_mm_repeat(X, wr, wc)
		return(X)
	}
	else {
		_mm_repeat(Xnew=X, wr, wc)
		return(Xnew)
	}
}

void _mm_repeat(transmorphic matrix X, real scalar wr,
 real scalar wc)
{
	real scalar i, rr, rc, r, c

	r = rows(X)
	c = cols(X)
	rr = max((trunc(wr-1),0))
	if (rr>0) {
		X = X \ J(rr*r, c, missingof(X))
		if (r>0 & c>0)
		 for (i=1; i<=rr; i++) X[|i*r+1,1 \ i*r+r,.|] = X[|1,1 \ r,.|]
	}
	rc = max((trunc(wc-1),0))
	if (rc>0) {
		X = X , J(rows(X), rc*c, missingof(X))
		if (r>0 & c>0)
		 for (i=1; i<=rc; i++) X[|1,i*c+1 \ .,i*c+c|] = X[|1,1 \ .,c|]
	}
}

end
