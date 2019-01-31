*! version 1.0.0, Ben Jann, 04jul2006
version 9.2

local opts "o1, o2, o3, o4, o5, o6, o7, o8, o9, o10"

mata:

struct mm_jkstats {
	real colvector reps, ndrop, fpc
	real rowvector stat
	real matrix    rstat
}

// note: reps is total n. of replications
//       ndrop is n. of unsuccessful replications
//       rows(rstat) = reps-ndrop is n. of successful replications

struct mm_jkstats scalar mm_jk(
 pointer(real rowvector function) scalar f,
 real matrix X,
 | real colvector w,
   real scalar nodots,       // 0
   real colvector strata,    // .
   real colvector cluster,   // .
   real colvector subpop,    // .
   real colvector fpc,       // .
   real matrix stat,         // .
   `opts')                   // opts to pass to f
{
	real scalar               i, ii, i0, i1, h, h0, h1, j
	real colvector            C, wjk, wjki
	real matrix               S
	transmorphic              fsetup
	struct mm_jkstats scalar  jk
	pointer scalar            ws

	if (args()<3) w = 1
	if (args()<4) nodots = 0
	if (args()<7) subpop = .
	if (args()<8) fpc = .
	if (args()<9) stat = .

	if (subpop!=. & rows(subpop)>0) ws = &(w :* (subpop:!=0))
	else ws = &w
	if (rows(*ws)==1) ws = &J(rows(X), 1, *ws)

	fsetup = mm_callf_setup(f, args()-9, `opts')

	if (stat==.) jk.stat = mm_callf(fsetup, X, *ws)[1,]
	else         jk.stat = stat[1,]

	mm_panels(strata, S=J(1,2,rows(X)), cluster, C=.)
	if (anyof(S, 1)) _error(460, "singleton cluster detected")

	if (nodots==0) {
		printf("{txt}\nJackknife replications ({res}%g{txt})\n",
		 colsum(S[,2]))
		display("{txt}{hline 4}{c +}{hline 3} 1 " +
		 "{hline 3}{c +}{hline 3} 2 " + "{hline 3}{c +}{hline 3} 3 " +
		 "{hline 3}{c +}{hline 3} 4 " + "{hline 3}{c +}{hline 3} 5 ")
	}
	jk.rstat = J(colsum(S[,2]), cols(jk.stat), .)
	jk.ndrop = J(rows(S), 1, 0)
	jk.fpc   = J(rows(S), 1, 0)
	h1 = i1 = 0
	i = ii = 1
	for (h=1; h<=rows(S); h++) { // for each stratum
		h0 = h1 + 1
		h1 = h1 + S[h,1]
		if (fpc!=. & rows(fpc)>0) jk.fpc[h] = fpc[h1]
		wjk = *ws
		wjk[| h0 \ h1 |] = wjk[| h0 \ h1 |] * (S[h,2] / (S[h,2] - 1))
		for (j=1; j<=S[h,2]; j++) { // for each PSU
			i0 = i1 + 1
			i1 = i1 + (C!=. ? C[i] : 1)
			wjki = wjk
			wjki[| i0 \ i1 |] = J(i1-i0+1,1,0)
			jk.rstat[ii,] = mm_callf(fsetup, X, wjki)[1,]
			if (missing(jk.rstat[ii,])) {
				jk.ndrop[h] = jk.ndrop[h] + 1
				if (nodots==0) printf("{err}x{txt}")
			}
			else {
				ii++
				if (nodots==0) printf(".")
			}
			if (nodots==0) {
				if (mod(i,50)==0) printf(" %5.0f\n",i)
				displayflush()
			}
			i++
		}
	}
	if (nodots==0 & mod(i-1,50)) display("")
	if (ii<i) {
		if (ii>1) jk.rstat = jk.rstat[|1,1 \ ii-1,.|]
		else jk.rstat = J(0, cols(jk.rstat), .)
	}
	jk.reps = S[,2]
	return(jk)
}

real matrix mm_jk_report(
 struct mm_jkstats scalar jk,    // jackknife estimates
 | string vector what,           // b,pseudo,mean,bias,v,se,ci
   real scalar level,            //
   real scalar mse,              // !=0 => use mse formula
   real vector fpc0)             // sampling fraction (for each stratum)
{
	real scalar    theta, pseudo, mean, bias, se, v, ci, l,
	               i, alpha, tval
	real rowvector jkse, jkv, jkmean
	real colvector fpc, nh, mh
	real matrix    res, jkpseudo

	if (args()<2) what  = "se"
	if (args()<3) level = st_numscalar("c(level)")
	if (args()<4) mse   = 0

	theta = pseudo = mean = bias = v = se = ci = l = 0
	for (i=1; i<=length(what); i++) {
		if      (what[i]=="theta" | what[i]=="b") theta = (l = l + 1)
		else if (what[i]=="pseudo") pseudo = (l = l + rows(jk.rstat))
		else if (what[i]=="mean"  ) mean  = (l = l + 1)
		else if (what[i]=="bias"  ) bias  = (l = l + 1)
		else if (what[i]=="v"     ) v     = (l = l + cols(jk.rstat))
		else if (what[i]=="se"    ) se    = (l = l + 1)
		else if (what[i]=="ci"    ) ci    = (l = l + 2)
		else _error(3498, "'" + what[i] + "' not allowed")
	}

	res = J(l, cols(jk.rstat), .)
	nh = jk.reps - jk.ndrop
	fpc = fpc0
	if (cols(fpc)!=1) fpc = fpc'
	if (rows(fpc)==0 | fpc==.) fpc = jk.fpc
	if (rows(fpc)==0)          fpc = J(rows(jk.reps),1,0)
// observed parameter vector
	if (theta) res[theta,] = jk.stat
// pseudo values
	if (pseudo | (mse==0 & (mean | bias))) {
		jkpseudo = jk.rstat +
		 mm_expand(nh, nh, 1, 1):*(jk.stat :- jk.rstat)
		if (pseudo)
		 res[|pseudo-rows(jk.rstat)+1,1 \ pseudo,.|] = jkpseudo
	}
// jackknife mean
	if (mean | bias) {
		if (mse) jkmean = mean(jk.rstat, 1)
		else     jkmean = mean(jkpseudo, 1)
		if (mean) res[mean,] = jkmean
	}
// bias
	if (bias) res[bias,] = jkmean - jk.stat
// variance
	if (v) {
		mh = mm_expand((1:-fpc) :* (nh :- 1) :/ nh, nh, 1, 1)
		if (mse) jkv = mm_sse(jk.rstat, mh, jk.stat)
		else     jkv = mm_sse(jk.rstat, mh, mean(jk.rstat,1))
		res[|v-cols(res)+1,1 \ v,.|] = jkv
	}
// standard errors
	if (se | ci) {
		if (v) jkse = sqrt(diagonal(jkv))'
		else {
			mh = mm_expand((1:-fpc) :* (nh :- 1) :/ nh, nh, 1, 1)
			if (mse) jkse = sqrt(mm_colsse(jk.rstat, mh, jk.stat))
			else     jkse = sqrt(mm_colsse(jk.rstat, mh, mean(jk.rstat,1)))
		}
		if (se) res[se,] = jkse
	}
// confidence interval
	if (ci) {
		alpha = (100-level)/200
		tval  = invttail(rows(jk.rstat)-rows(jk.reps), alpha)
		res[ci-1,] = jk.stat - tval*jkse
		res[ci,]   = jk.stat + tval*jkse
	}
	return(res)
}

end
