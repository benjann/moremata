*! version 1.0.0, Ben Jann, 10jul2006
version 9.2

local opts "o1, o2, o3, o4, o5, o6, o7, o8, o9, o10"

mata:

struct mm_bsstats {
	real scalar    reps, ndrop
	real rowvector stat, se
	real matrix    rstat, rse
}

// note: reps is total n. of replications
//       ndrop is n. of unsuccessful replications
//       rows(rstat) = reps-ndrop is n. of successful replications

struct mm_bsstats scalar mm_bs(
 pointer(real rowvector function) scalar f,
 real matrix X,
 | real colvector w,
   real scalar reps,        // 50
   real scalar d,           // 0
   real scalar nodots,      // 0
   real colvector strata,   // .
   real colvector cluster,  // .
   real matrix stat,        // .
   `opts')                  // opts to pass to f
{
	real scalar               i, ii, hasse
	real colvector            C, N, s
	real matrix               S, res
	transmorphic              setup
	struct mm_bsstats scalar  bs

	if (args()<3) w = 1
	if (args()<4) reps = 50
	if (args()<5) nodots = 0
	if (args()<6) d = 0
	if (args()<9) stat = .

	setup = mm_callf_setup(f, args()-9, `opts')

	if (stat==.) {
		res = mm_callf(setup, X, w)
		bs.stat = res[1,]
		if (hasse=(rows(res)>1)) bs.se = res[2,]
	}
	else {
		bs.stat = stat[1,]
		if (hasse=(rows(stat)>1)) bs.se = stat[2,]
	}

	mm_panels(strata, S=J(1,2,rows(X)), cluster, C=.)
	if (anyof(S, 1)) _error(460, "singleton cluster detected")
	if (d==0|d>=.) N = .
	else N = S[,2] :- d

	if (nodots==0) {
		printf("{txt}\nBootstrap replications ({res}%g{txt})\n", reps)
		display("{txt}{hline 4}{c +}{hline 3} 1 " +
		 "{hline 3}{c +}{hline 3} 2 " + "{hline 3}{c +}{hline 3} 3 " +
		 "{hline 3}{c +}{hline 3} 4 " + "{hline 3}{c +}{hline 3} 5 ")
	}
	bs.rstat = J(reps, cols(bs.stat), .)
	if (hasse) bs.rse = J(reps, cols(bs.se), .)
	for (i=ii=1; i<=reps; i++) {
		s  = mm_sample(N, S, C)
		res = mm_callf(setup, X[s,], rows(w)==1 ? w : w[s,])
		if (missing(res)) {
			if (nodots==0) printf("{err}x{txt}")
		}
		else {
			bs.rstat[ii,] = res[1,]
			if (hasse) bs.rse[ii,] = res[2,]
			ii++
			if (nodots==0) printf(".")
		}
		if (nodots==0) {
			if (mod(i,50)==0) printf(" %5.0f\n",i)
			displayflush()
		}
	}
	if (nodots==0 & mod(i-1,50)) display("")
	if (ii<i) {
		if (ii>1) {
			bs.rstat = bs.rstat[|1,1 \ ii-1,.|]
			if (hasse) bs.rse = bs.rse[|1,1 \ ii-1,.|]
		}
		else {
			bs.rstat = J(0, cols(bs.rstat), .)
			if (hasse) bs.rse = J(0, cols(bs.rse), .)
		}
	}
	bs.reps = reps
	bs.ndrop = reps - rows(bs.rstat)
	return(bs)
}

struct mm_bsstats scalar mm_bs2(
 pointer(real rowvector function) scalar f,
 real matrix X,
 | real colvector w,
   real scalar reps,        // 50
   real scalar d,           // 0
   real scalar nodots,      // 0
   real colvector strata,   // .
   real colvector cluster,  // .
   real matrix stat,        // .
   `opts')                  // opts to pass to f
{
	real scalar               i, ii, hasse
	real colvector            C, N, ws
	real matrix               S, res
	transmorphic              setup
	struct mm_bsstats scalar  bs

	if (args()<3) w = 1
	if (args()<4) reps = 50
	if (args()<5) nodots = 0
	if (args()<6) d = 0
	if (args()<9) stat = .

	setup = mm_callf_setup(f, args()-9, `opts')

	if (stat==.) {
		res = mm_callf(setup, X, w)
		bs.stat = res[1,]
		if (hasse=(rows(res)>1)) bs.se = res[2,]
	}
	else {
		bs.stat = stat[1,]
		if (hasse=(rows(stat)>1)) bs.se = stat[2,]
	}

	mm_panels(strata, S=J(1,2,rows(X)), cluster, C=.)
	if (anyof(S, 1)) _error(460, "singleton cluster detected")
	if (d==0|d>=.) N = .
	else N = S[,2] :- d

	if (nodots==0) {
		printf("{txt}\nBootstrap replications ({res}%g{txt})\n", reps)
		display("{txt}{hline 4}{c +}{hline 3} 1 " +
		 "{hline 3}{c +}{hline 3} 2 " + "{hline 3}{c +}{hline 3} 3 " +
		 "{hline 3}{c +}{hline 3} 4 " + "{hline 3}{c +}{hline 3} 5 ")
	}
	bs.rstat = J(reps, cols(bs.stat), .)
	if (hasse) bs.rse = J(reps, cols(bs.se), .)
	for (i=ii=1; i<=reps; i++) {
		ws = mm_sample(N, S, C, 1, 0, 1):*w
		res = mm_callf(setup, X, ws)
		if (missing(res)) {
			if (nodots==0) printf("{err}x{txt}")
		}
		else {
			bs.rstat[ii,] = res[1,]
			if (hasse) bs.rse[ii,] = res[2,]
			ii++
			if (nodots==0) printf(".")
		}
		if (nodots==0) {
			if (mod(i,50)==0) printf(" %5.0f\n",i)
			displayflush()
		}
	}
	if (nodots==0 & mod(i-1,50)) display("")
	if (ii<i) {
		if (ii>1) {
			bs.rstat = bs.rstat[|1,1 \ ii-1,.|]
			if (hasse) bs.rse = bs.rse[|1,1 \ ii-1,.|]
		}
		else {
			bs.rstat = J(0, cols(bs.rstat), .)
			if (hasse) bs.rse = J(0, cols(bs.rse), .)
		}
	}
	bs.reps = reps
	bs.ndrop = reps - rows(bs.rstat)
	return(bs)
}

real matrix mm_bs_report(
 struct mm_bsstats scalar bs, // boostrap estimates
 | string vector what,        // b,mean,bias,v,se,n,basic,p,bc,bca,t
   real scalar level,         //
   real scalar mse,           // !=0 => use mse formula
   struct mm_jkstats jk)      // jackknife estimates
{
	real scalar    theta, mean, bias, v, se, n, basic, p, bc,
	               bca, t, l, i, alpha, tval
	real rowvector bsv, bsse, bsmean, z0, a
	real matrix    res

	if (args()<2) what  = "se"
	if (args()<3) level = st_numscalar("c(level)")
	if (args()<4) mse   = 0

	theta = mean = bias = v = se = n = basic = p = bc = bca = t = l = 0
	for (i=1; i<=length(what); i++) {
		if      (what[i]=="theta" | what[i]=="b") theta = (l = l + 1)
		else if (what[i]=="mean" ) mean  = (l = l + 1)
		else if (what[i]=="bias" ) bias  = (l = l + 1)
		else if (what[i]=="v"    ) v     = (l = l + cols(bs.rstat))
		else if (what[i]=="se"   ) se    = (l = l + 1)
		else if (what[i]=="n" | what[i]=="ci") n = (l = l + 2)
		else if (what[i]=="basic") basic = (l = l + 2)
		else if (what[i]=="p"    ) p     = (l = l + 2)
		else if (what[i]=="bc"   ) bc    = (l = l + 2)
		else if (what[i]=="bca"  ) bca   = (l = l + 2)
		else if (what[i]=="t"    ) t     = (l = l + 2)
		else _error(3498, "'" + what[i] + "' not allowed")
	}

	res = J(l, cols(bs.rstat), .)
// observed parameter vector
	if (theta) res[theta,] = bs.stat
// mean of bs replicates
	if (mean | bias) {
		bsmean = mean(bs.rstat, 1)
		if (mean) res[mean,] = bsmean
	}
// bias
	if (bias) res[bias,] = bsmean - bs.stat
// variance
	if (v) {
		if (mse) bsv = mm_mse(bs.rstat, 1, bs.stat)
		else     bsv = variance(bs.rstat, 1)
		res[|v-cols(res)+1,1 \ v,.|] = bsv
	}
// standard errors
	if (se | n) {
		if (v) bsse = sqrt(diagonal(bsv))'
		else {
			if (mse) bsse = sqrt(mm_colmse(bs.rstat, 1, bs.stat))
			else     bsse = sqrt(mm_colvar(bs.rstat, 1))
		}
		if (se) res[se,] = bsse
	}
// normal ci
	alpha = (100-level)/200
	if (n) {
		tval = invttail(bs.reps-bs.ndrop-1, alpha)
		res[n-1,] = bs.stat - tval*bsse
		res[n,]   = bs.stat + tval*bsse
	}
// basic ci
	if (basic) res[|basic-1,1 \ basic,.|] = 2*bs.stat :-
	 mm_quantile(bs.rstat, 1, ((1-alpha) \ alpha))
// percentile ci
	if (p) res[|p-1,1 \ p,.|] =
	 mm_quantile(bs.rstat, 1, (alpha \ (1-alpha)))
// bias-corrected ci
	if (bc | bca) z0 = editmissing(invnormal(
	 colsum(bs.rstat:<=bs.stat)/rows(bs.rstat) ), 0)
	if (bc) {
		res[|bc-1,1 \ bc,.|] = mm_quantile(bs.rstat, 1,
		 (normal(2*z0 :+ invnormal(alpha)) \
		 normal(2*z0 :- invnormal(alpha))))
	}
// bias-corrected and accelerated ci
	if (bca) {
		a = colsum((mean(jk.rstat, 1):-jk.rstat):^3) :/
		    ( 6 * colsum((mean(jk.rstat, 1):-jk.rstat):^2):^1.5 )
		res[|bca-1,1 \ bca,.|] = mm_quantile(bs.rstat, 1,
		 (normal(z0 :+ (z0 :+ invnormal(alpha)) :/
		   (1:-a:*(z0 :+ invnormal(alpha)))) \
		 normal(z0 :+ (z0 :- invnormal(alpha)) :/
		   (1:-a:*(z0 :- invnormal(alpha))))))
	}
// percentile-t ci
	if (t) {
		res[|t-1,1 \ t,.|] = bs.stat :- bs.se :* mm_quantile(
		 editmissing((bs.rstat:-bs.stat):/bs.rse, 0), 1,
		 ((1-alpha) \ alpha) )
	}
	return(res)
}

end
