*! version 1.0.1, Ben Jann, 14apr2006
version 9.1
mata:

real colvector mm_sample(
 real colvector n,
 real matrix strata,
 | real colvector cluster,
   real colvector w,
   real scalar wor,
   real scalar count,
   real scalar fast)
{
	real scalar i, b, e, n0, n1, N, cb, ce
	real colvector s, si, ci, wi, nn, R
	pointer scalar f

	if (args()<3) cluster=.
	if (args()<4) w=1
	if (args()<5) wor=0
	if (args()<6) count=0
	if (args()<7) fast=0
	if (cols(strata)<1) _error(3200)
	if (fast==0) {
		if (rows(strata)>1 | (cluster==. & rows(w)==1)) {
			if (missing(strata)) _error(3351, "'strata' has missing values")
		}
	}
	if (rows(n)<1) _error(3200)

//choose sampling function and check w
	if (rows(w)==1) {
		if (wor) f = &mm_srswor()
		else     f = &mm_srswr()
	}
	else {
		if (cluster!=. & rows(w)!=rows(cluster))
		  _error(3200, "rows('w') unequal number of clusters")
		if (fast==0) {
			if (missing(w)) _error(3351, "'w' has missing values")
			if (colsum(w:<0)) _error(3498, "'w' has negative values")
			if (colsum(w)<=0) _error(3498, "sum('w') is zero")
		}
		if (wor) f = &mm_upswor()
		else     f = &mm_upswr()
	}

//simple sample, unstratified
	if (rows(strata)==1) {
		N = strata[1,1]
		if (cluster==.) {
			if (N>=.) N = rows(w)
			else if (rows(w)!=1 & N!=rows(w))
			  _error(3200, "rows('w') unequal population size")
			return((*f)(n, (rows(w)==1 ? N : w), count))
		}

//cluster sample, unstratified
		if (rows(cluster)<1) _error(3200)
		if (N>=.) N = colsum(cluster)
		else if (fast==0) {
			if (N!=colsum(cluster))
			  _error(3200, "sum of cluster sizes unequal population size")
		}
		s = (*f)(n, (rows(w)==1 ? rows(cluster) : w), count)
		return(_mm_expandclusters(s, cluster, count, N))
	}

//stratified sample
	if (rows(strata)<1) _error(3200)
	if (cluster!=. & cols(strata)<2) _error(3200)
//-sample sizes
	if (rows(n)==1) {
		if (n>=.) {
			if (cluster==.) nn = strata[.,1]
			else            nn = strata[.,2]
		}
		else nn = J(rows(strata),1,n)
	}
	else {
		if (rows(n)!=rows(strata)) _error(3200)
		nn = n
		for (i=1;i<=rows(nn);i++) {
			if (nn[i]>=.) {
				if (cluster==.) nn[i] = strata[i,1]
				else            nn[i] = strata[i,2]
			}
		}
	}
//-prepare sample vector
	if (count) s = J(colsum(strata[.,1]),1,.)
	else s = J(colsum(nn),1,.)
//-draw samples within strata
	if (count==0) n0 = 1
	e = 0
	ce = 0
	for (i=1;i<=rows(strata);i++) {
		b = e + 1
		e = e + strata[i,1]
//--simple sample
		if (cluster==.) {
			si = (*f)(nn[i], (rows(w)==1 ? strata[i,1] : w[|b \ e|]), count)
		}
//--cluster sample
		else {
			cb = ce + 1
			ce = ce + strata[i,2]
			ci = cluster[|cb \ ce|]
			wi = (rows(w)==1 ? rows(ci) : w[|cb \ ce|])
			if (fast==0) {
				if (strata[i,1]!=colsum(ci))
				  _error(3200, "sum of cluster sizes unequal size of stratum")
			}
			si = (*f)(nn[i], wi, count)
			if (count) si = _mm_expandclusters(si, ci, 1, strata[i,1])
		}
//--add subsample to sample vector
		if (count) R = (b \ e)
		else {
			n1 = n0 + rows(si) - 1
			R = (n0 \ n1)
			if (cluster==.) si = si :+ (b-1)
			else si = si :+ (cb-1)
			n0 = n1 + 1
		}
		if (R[1]<=R[2]) s[|R|] = si
	}
	if (count==0&cluster!=.) s = _mm_expandclusters(s, cluster)
	return(s)
}

real colvector _mm_expandclusters(real colvector s0,
 real colvector cluster, | real scalar count, real scalar N)
{
	real scalar i, e, b, n0, n1
	real colvector s, eclust

	if (args()<3) count = 0
	if (count) {
		s = J(N,1,.)
		e = 0
		for (i=1;i<=rows(cluster);i++) {
			b = e + 1
			e = e + cluster[i]
			s[|b \ e|] = J(e-b+1, 1, s0[i])
		}
		return(s)
	}
	s = J(colsum(cluster[s0,]),1,.)
	eclust = mm_colrunsum(cluster)
	n0 = 1
	for (i=1;i<=rows(s0);i++) {
		e = eclust[s0[i]]
		b = e - cluster[s0[i]] + 1
		n1 = n0 + e-b
		s[|n0 \ n1|] = (b::e)
		n0 = n1+1
	}
	return(s)
}

end
