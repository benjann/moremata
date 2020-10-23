*! version 1.0.2, Ben Jann, 23oct2020
version 9.1
mata:

real colvector mm_sample(
 real colvector n,
 real matrix strata,
 | real colvector cluster,
   real colvector w,
   real scalar wor,
   real scalar count,
   real scalar fast,
   real scalar alt,
   real scalar nowarn)
{
	real scalar i, b, e, n0, n1, N, cb, ce, ups
	real colvector s, si, ci, wi, nn, R

	if (args()<3) cluster=.
	if (args()<4) w=1
	if (args()<5) wor=0
	if (args()<6) count=0
	if (args()<7) fast=0
	if (args()<8) alt=0
	if (args()<9) nowarn=0
	if (cols(strata)<1) _error(3200)
	ups = (rows(w)!=1)
	if (fast==0) {
		if (rows(strata)>1 | (cluster==. & ups==0)) {
			if (missing(strata)) _error(3351, "'strata' has missing values")
		}
	}
	if (rows(n)<1) _error(3200)

//check w
	if (ups) {
		if (cluster!=.) {
			if (rows(w)!=rows(cluster))
				_error(3200, "rows('w') unequal number of clusters")
		}
		if (fast==0) {
			if (missing(w)) _error(3351, "'w' has missing values")
			if (colsum(w:<0)) _error(3498, "'w' has negative values")
			if (colsum(w)<=0) _error(3498, "sum('w') is zero")
		}
	}

//simple sample, unstratified
	if (rows(strata)==1) {
		N = strata[1,1]
		if (cluster==.) {
			if (N>=.) N = rows(w)
			else if (ups) {
				if (N!=rows(w))
					_error(3200, "rows('w') unequal population size")
			}
			return(_mm_sample(n, (ups ? w : N), ups, wor, count, alt, nowarn))
		}

//cluster sample, unstratified
		if (rows(cluster)<1) _error(3200)
		if (N>=.) N = colsum(cluster)
		else if (fast==0) {
			if (N!=colsum(cluster))
			  _error(3200, "sum of cluster sizes unequal population size")
		}
		s = _mm_sample(n, (ups ? w : rows(cluster)), ups, wor, count, alt, nowarn)
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
			si = _mm_sample(nn[i], (ups ? w[|b \ e|] : strata[i,1]), ups, wor,
				count, alt, nowarn)
		}
//--cluster sample
		else {
			cb = ce + 1
			ce = ce + strata[i,2]
			ci = cluster[|cb \ ce|]
			wi = (ups ? w[|cb \ ce|] : rows(ci))
			if (fast==0) {
				if (strata[i,1]!=colsum(ci))
				  _error(3200, "sum of cluster sizes unequal size of stratum")
			}
			si = _mm_sample(nn[i], wi, ups, wor, count, alt, nowarn)
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

real colvector _mm_sample(real scalar n, real colvector NorW, 
    real scalar ups, real scalar wor, real scalar count, real scalar alt, 
    real scalar nowarn)
{
    if (ups) {
        if (wor) return(mm_upswor(n, NorW, count, nowarn))
        return(mm_upswr(n, NorW, count))
    }
    if (wor) return(mm_srswor(n, NorW, count, alt))
    return(mm_srswr(n, NorW, count))
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
