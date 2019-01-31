*! version 1.0.1, Ben Jann, 02may2007
version 9.2
mata:

/*-SUBSETS---------------------------------------------------------------*/
// Generate all combinations (subsets) one-by-one (lexicographic)
// based on Algorithm 5.8 from Reingold et al. (1977)
// Usage:
//
//        info = mm_subsetsetup(n,k)
//        while ((c = mm_subset(info)) != J(0,1,.)) {
//             ... c ...
//        }
//
// The total number of combinations can be computed using
// Mata's -comb()- function, i.e.
//
//        Ctot = comb(n, k)

struct mm_subsetinfo {
	real scalar    n, k, i, j, counter, algorithm
	real colvector x, y
}

struct mm_subsetinfo scalar mm_subsetsetup(
 real scalar n,
 | real scalar k)
{
	struct mm_subsetinfo scalar s
// setup
	s.n = trunc(n)
	s.k = trunc(k)
	if (s.n>=.) _error(3351)
	if (s.k>=.) s.k = s.n
	else if (s.k>s.n) _error(3300,"k may not be larger than n")
	if (s.n<1|s.k<1) s.i = 0  // => return empty set
	else {
		s.x = -1 \ (1::s.k)
		s.i = 1
	}
	s.i = s.i + 1 // offset i by 1
	s.k = s.k + 1 // offset k by 1
	s.counter = 0
	return(s)
}

real colvector mm_subset(struct mm_subsetinfo scalar s)
{
	real scalar    j
	real colvector subset

	if (s.i==1) return(J(0,1,.)) /*done*/
// get current subset (to be returned at end)
	subset = s.x[|2 \ s.k|]
// generate next subset
	s.i = s.k
	while (s.x[s.i] == s.n-s.k+s.i) {
		s.i = s.i - 1
	}
	s.x[s.i] = s.x[s.i] + 1
	for (j=s.i+1;j<=s.k;j++) s.x[j] = s.x[j-1] + 1
/* the following would improve speed but does not work for some reason
	if (s.i<s.k) s.x[|s.i+1 \ s.k|] = s.x[|s.i \ s.k-1|] :+ 1
*/
	s.counter = s.counter + 1
	return(subset)
}

// wrapper for mm_subsetinfo()/mm_subset()
real matrix mm_subsets(real scalar n, | real scalar k)
{
	real scalar    i
	real colvector set
	real matrix    res
	struct mm_subsetinfo scalar s

	s = mm_subsetsetup(n, (k<. ? k : n))
	if ((set = mm_subset(s)) == J(0,1,.)) return(set)
	res = J(rows(set), comb(n, rows(set)), .)
	res[,i=1] = set
	while ((set = mm_subset(s)) != J(0,1,.)) res[,++i] = set
	return(res)
}


/*-COMPOSITIONS----------------------------------------------------------*/
// Generate all k-way compositions of n, one-by-one
// Default algorithm: direct (less elegant but faster) (anti-lexicographic)
// Alternate algorithm: indirect via subsets (lexicographic)
//
// Usage:
//
//        info = mm_compositionsetup(n,k)
//        while ((c = mm_composition(info)) != J(0,1,.)) {
//             ... c ...
//        }
//
// or:
//
//        info = mm_compositionsetup(n,k,1)
//        while ((c = mm_composition(info)) != J(0,1,.)) {
//             ... c ...
//        }
//
// mm_ncompositions() returns the total number of compositions

real scalar mm_ncompositions(
 real scalar n,
 | real scalar k0)
{
	real scalar k

	k = trunc(k0<. ? k0 :  n)
	return(comb(trunc(n)+k-1,k-1))
}

struct mm_subsetinfo scalar mm_compositionsetup(
 real scalar n,
 | real scalar k,
   real scalar alt) // !=0 => alternate algorithm
{
	struct mm_subsetinfo scalar s

	if (args()<3) alt = 0
	if (alt) {
		s = _mm_composition2setup(n,k)
		s.algorithm = 2
		return(s)
	}
	s = _mm_compositionsetup(n,k)
	s.algorithm = 1
	return(s)
}

real colvector mm_composition(struct mm_subsetinfo scalar s)
{
	if (s.algorithm==1) return(_mm_composition(s))
	if (s.algorithm==2) return(_mm_composition2(s))
	_error(3498, "unknown algorithm")
}

// default algorithm

struct mm_subsetinfo scalar _mm_compositionsetup(
 real scalar n,
 | real scalar k)
{
	struct mm_subsetinfo scalar s

	s.n = trunc(n)
	s.k = (k<. ? trunc(k) :  s.n)
	if (s.n>=.) _error(3351)
	if (s.k<1) {
		s.i = 0
		s.x = J(0,1,.)
	}
	else if (s.n<1) {
		s.i = 0
		s.x = J(s.k,1,0)
	}
	else if (s.k<2) {
		s.i = 0
		s.x = J(s.k,1,s.n)
	}
	else {
		s.x = s.n \ J(s.k-1,1,0)
		s.i = 1
	}
	s.j = 2
	s.counter = 0
	return(s)
}

real colvector _mm_composition(struct mm_subsetinfo scalar s)
{
	real colvector subset

	subset = s.x
	if (s.i<1) {
		if (subset!=J(0,1,.)) s.counter = s.counter + 1
		s.x = J(0,1,.)
		return(subset) /*done*/
	}
	while (s.x[s.i]==0 & s.i>1) { // move back
		s.x[s.i] = s.x[s.j]
		s.x[s.j] = 0
		s.j = s.i--
	}
	if (s.x[s.i]>0) {
		s.x[s.i] = s.x[s.i]-1
		s.x[s.j] = s.x[s.j]+1
		if (s.j<s.k) s.i = s.j++
	}
	else {
		s.i = 0
		s.x = J(0,1,.)
	}
	s.counter = s.counter + 1
	return(subset)
}

// alternate algorithm

struct mm_subsetinfo scalar _mm_composition2setup(
 real scalar n,
 | real scalar k0)
{
	real scalar k

	k = trunc(k0<. ? k0 :  n)
	return(mm_subsetsetup(trunc(n)+k-1,k-1))
}

real colvector _mm_composition2(struct mm_subsetinfo scalar s)
{
	real colvector subset

	if (s.i==1) {
		if (s.k==1) {  // length of composition is 1
			s.k = 0
			s.counter = s.counter + 1
			return(s.n)
		}
		return(J(0,1,.)) /*done*/
	}
	subset = mm_subset(s)
	return( (subset \ s.n+1) - (0 \ subset) :- 1 )
}

// wrapper for mm_compositioninfo()/mm_composition()
real matrix mm_compositions(real scalar n, | real scalar k, real scalar alt)
{
	real scalar    i
	real colvector set
	real matrix    res
	struct mm_subsetinfo scalar s

	if (args()<3) alt = 0
	s = mm_compositionsetup(n, (k<. ? k : n), alt)
	if ((set = mm_composition(s)) == J(0,1,.)) return(set)
	res = J(rows(set), mm_ncompositions(n, rows(set)), .)
	res[,i=1] = set
	while ((set = mm_composition(s)) != J(0,1,.)) res[,++i] = set
	return(res)
}

/*-RANDOM SUBSETS--------------------------------------------------------*/
// Return random combination (subsets)
// Usage:
//
//        for (i=1;i<=reps;i++) {
//             c = mm_rsubset(n, k)
//             ... c ...
//        }
//
// mm_rsubset() returns jumbled sets; apply sort() if you need
// ordered sets, e.g.
//
//        c = sort(mm_rsubset(n, k), 1)

real colvector mm_rsubset(
 real scalar n,
 | real scalar k)
{
	if (k<. & k>n) _error(3300, "k may not be larger than n")
	if (n<1|k<1) return(J(0,1,.))
	return(mm_unorder2(n)[|1 \ (k<. ? k : n)|])
}

/*-RANDOM COMPOSITIONS---------------------------------------------------*/
// Return random combination (subsets) and random composition
// Usage:
//
//        for (i=1;i<=reps;i++) {
//             c = mm_rcomposition(n, k)
//             ... c ...
//        }

real colvector mm_rcomposition(
 real scalar n0,
 | real scalar k0)
{
	real scalar    n, k
	real colvector subset

	n = trunc(n0)
	k = (k0<. ? trunc(k0) : n)
	if (n>=.) _error(3351)
	if (k<1) return(J(0,1,.))
	if (n<1) return(J(k,1,0))
	if (k<2) return(J(k,1,n))
	subset = sort(mm_rsubset(n+k-1,k-1), 1)
	return( (subset \ n+k) - (0 \ subset) :- 1 )
}

/*-PARTITIONS------------------------------------------------------------*/
// Generate all k-way partitions of n, one-by-one
// Default algorithm: based on ZS1 algorithm in Zoghbi&Stojmenovic (1998)
//                    (anti-lexicographic)
// Alternate algorithm: based on Algorithm 5.12 in Reingold et al. (1977)
//                    (dictionary order)
// Usage:
//
//        info = mm_partitionsetup(n,k)
//        while ((c = mm_partition(info)) != J(0,1,.)) {
//             ... c ...
//        }
//
// or:
//
//        info = mm_partitionsetup(n,k,1)
//        while ((c = mm_partition(info)) != J(0,1,.)) {
//             ... c ...
//        }
//
// mm_npartitions() returns the total number of compositions, based
// on algorithms from
// http://home.att.net/~numericana/answer/numbers.htm#partitions

// number of partitions of n with k or fewer addends (slow)
real scalar mm_npartitions(
 real scalar n,
 | real scalar k)
{
	real scalar      u, i, j, m
	real colvector   p, a

	m = min((n,k))
	if (m==n) return(_mm_npartitions(n))
	p = J(n+1,1,1)
	for (u=2;u<=m;u++) {
		a = p
		p = J(n+1,1,0)
		for (i=0;i<=n;i=i+u) {
			for (j=i+1;j<=n+1;j++) {
				p[j] = p[j] + a[j-i]
			}
		}
	}
	return(p[n+1])
}
// total number of partitions of n (fast)
real scalar _mm_npartitions(
 real scalar n)
{
	real scalar      i, j, s, k
	real colvector   p

	p = J(n+1,1,1)
	for (i=1;i<=n;i++) {
		j = 1 ; k = 1 ; s = 0
		while (j>0) {
			j = i - (3*k*k + k)/2
			if (j>=0) s = s - (-1)^k * p[j+1]
			j = i - (3*k*k - k)/2
			if (j>=0) s = s - (-1)^k * p[j+1]
			k = k + 1
		}
		p[i+1] = s
	}
	return(p[n+1])
}

struct mm_subsetinfo scalar mm_partitionsetup(
 real scalar n,
 | real scalar k,
   real scalar pad,  // zero-padded output
   real scalar alt)  // !=0 => alternate algorithm
{
	struct mm_subsetinfo scalar s

	if (args()<3) pad = 0
	if (args()<4) alt = 0
	if (alt) {
		s = _mm_partition2setup(n,k)
		s.algorithm = 3 + (pad!=0)
		return(s)
	}
	s = _mm_partitionsetup(n,k)
	s.algorithm = 1 + (pad!=0)
	return(s)
}

real colvector mm_partition(struct mm_subsetinfo scalar s)
{
	if (s.algorithm==1) return(_mm_partition(s,0))
	if (s.algorithm==2) return(_mm_partition(s,1))
	if (s.algorithm==3) return(_mm_partition2(s,0))
	if (s.algorithm==4) return(_mm_partition2(s,1))
	_error(3498, "unknown algorithm")
}

// default algorithm

struct mm_subsetinfo scalar _mm_partitionsetup(
 real scalar n,
 | real scalar k)
{
	struct mm_subsetinfo scalar s

	s.n = trunc(n)
	s.k = trunc(k)
	if (s.n>=.) _error(3351)
	if (s.k>=.) s.k = s.n
//	else if (s.k>s.n) _error(3300,"k may not be larger than n")
	if (s.n<1 | s.k<1 ) s.i = -1
	else {
		s.x = s.n \ J(s.k-1, 1, 1)
		s.i = 1
	}
	s.j = 1
	s.counter = 0
	return(s)
}

real colvector _mm_partition(        // original ZS1 Algorithm does not
 struct mm_subsetinfo scalar s,      // support k-way partitions; changes
 real scalar pad)                    // are indicated
{
	real scalar    r, t, l
	real colvector subset

	if (s.i==-1) return(J(0,1,.)) /*done*/
	subset = s.x[|1 \ s.j|]
	if (pad) subset = subset \ J(s.k-rows(subset),1,0)
	l = (s.k<s.n ? s.k : s.n)     // added
	if (s.j>=l) {                 // added; original stopping rule is (s.x[1]==1)
		if (s.x[1]<=s.x[l]+1) {      // changed
			s.i = -1
			s.counter = s.counter + 1
			return(subset)
		}
		while (s.x[s.i]<=s.x[l]+1) {  // step back loop added
			s.j = s.j + s.x[s.i] - 1
			s.i = s.i - 1
		}
	}
	if (s.x[s.i]==2) {
		s.j      = s.j + 1
		s.x[s.i] = 1
		s.x[s.j] = 1             // added
		s.i      = s.i - 1
	}
	else {
		r = s.x[s.i] - 1
		t = s.j - s.i + 1
		s.x[s.i] = r
		while (t>=r) {
			s.i = s.i + 1
			s.x[s.i] = r
			t = t - r
		}
		if (t==0) {
			s.j = s.i
		}
		else {
			s.j = s.i + 1
			if (t==1) s.x[s.i+1] = 1    // added
			if (t>1) {
				s.i    = s.i + 1
				s.x[s.i] = t
			}
		}
	}
	s.counter = s.counter + 1
	return(subset)
}

// alternate algorithm

struct mm_subsetinfo scalar _mm_partition2setup(
 real scalar n,
 | real scalar k)
{
	real scalar offset
	struct mm_subsetinfo scalar s

	s.n = trunc(n)
	s.k = trunc(k)
	if (s.n>=.) _error(3351)
	if (s.k>=.) s.k = s.n
	offset = 2
	if (s.n<1 | s.k<1 ) s.i = 0 + offset
	else {
		s.i = 1 + offset
		s.x = s.y = J(s.n+offset, 1, .)
		s.x[-1+offset] = 0
		s.y[-1+offset] = 0
		s.x[0+offset]  = 0
		s.y[0+offset]  = s.n + 1
		s.x[1+offset]  = s.n
		s.y[1+offset]  = 1
	}
	s.counter = 0
	return(s)
}

real colvector _mm_partition2(
 struct mm_subsetinfo scalar s,
 real scalar pad)
{
	real scalar    skip, sum
	real colvector subset

	skip = 1
	while (skip) {
		if (s.i<=2) return(J(0,1,.)) /*done*/
		if ((skip = (sum(s.x[|3 \ s.i|])>s.k))==0)
		 subset = _mm_partition2_expand(s.x[|3 \ s.i|], s.y[|3 \ s.i|],
		 (pad ? s.k : .))
		sum = s.x[s.i]*s.y[s.i]
		if (s.x[s.i]==1) {
			s.i = s.i - 1
			sum = sum + s.x[s.i]*s.y[s.i]
		}
		if (s.y[s.i-1]==s.y[s.i]+1) {
			s.i      = s.i - 1
			s.x[s.i] = s.x[s.i] + 1
		}
		else {
			s.y[s.i] = s.y[s.i] + 1
			s.x[s.i] = 1
		}
		if (sum>s.y[s.i]) {
			s.y[s.i+1] = 1
			s.x[s.i+1] = sum - s.y[s.i]
			s.i        = s.i + 1
		}
	}
	s.counter = s.counter + 1
	return(subset)
}

real colvector _mm_partition2_expand(
 real vector m,
 real vector p,
 real scalar k)
{
	real scalar    i, j, l
	real colvector res

	if (k<.) res = J(k,1,0)
	else     res = J(sum(m),1,0)
	l = 0
	for (i=1;i<=rows(p);i++) {
		for (j=1;j<=m[i];j++) {
			res[++l] = p[i]
		}
	}
	return(res)
}

// wrapper for mm_partitioninfo()/mm_partition()
real matrix mm_partitions(real scalar n, | real scalar k, real scalar alt)
{
	real scalar    i
	real colvector set
	real matrix    res
	struct mm_subsetinfo scalar s

	if (args()<3) alt = 0
	s = mm_partitionsetup(n, (k<. ? k : n), 1, alt)
	if ((set = mm_partition(s)) == J(0,1,.)) return(set)
	res = J(rows(set), mm_npartitions(n, rows(set)), .)
	res[,i=1] = set
	while ((set = mm_partition(s)) != J(0,1,.)) res[,++i] = set
	return(res)
}

end
