# moremata
Stata module providing various Mata functions

  * `mm_kern()`: various kernel functions
  * `mm_kint()`: kernel integral functions
  * `mm_kderiv()`: kernel derivative functions
  * `mm_kdel0()`: canonical bandwidth of kernel
  * `mm_quantile()`: compute quantiles
  * `mm_median()`: compute median
  * `mm_iqrange()`: compute inter-quartile range
  * `mm_ecdf()`: compute cumulative distribution function
  * `mm_ecdf2()`: cumulative distribution at unique values
  * `mm_ranks()`: compute ranks/cumulative frequencies
  * `mm_relrank()`: compute relative ranks (grade transformation)
  * `mm_density()`: compute density (Stata 11 required)
  * `mm_ddens()`: compute density by diffusion
  * `mm_freq()`: compute frequency counts
  * `mm_histogram()`: produce histogram data
  * `mm_mgof()`: multinomial goodness-of-fit tests
  * `mm_collapse()`: summary statistics by subgroups
  * `_mm_collapse()`: summary statistics by subgroups
  * `mm_gini()`: Gini coefficient
  * `mm_nobs()`: number of observations
  * `mm_sample()`: draw random sample
  * `mm_srswr()`: SRS with replacement
  * `mm_srswor()`: SRS without replacement
  * `mm_upswr()`: UPS with replacement
  * `mm_upswor()`: UPS without replacement
  * `mm_bs()`: bootstrap estimation
  * `mm_bs2()`: bootstrap estimation
  * `mm_bs_report()`: report bootstrap results
  * `mm_jk()`: jackknife estimation
  * `mm_jk_report()`: report jackknife results
  * `mm_subset()`: obtain subsets, one at a time
  * `mm_composition()`: obtain compositions, one by one
  * `mm_ncompositions()`: determine number of compositions
  * `mm_partition()`: obtain partitions, one at a time
  * `mm_npartitionss()`: determine number of partitions
  * `mm_rsubset()`: draw random subset
  * `mm_rcomposition()`: draw random composition
  * `mm_greedy()`: one-to-one and one-to-many matching w/o replacement
  * `mm_greedy2()`: like `mm_greedy()`, but returning edge-list
  * `mm_greedy_pairs()`: transform result from `mm_greedy()` into edge-list
  * `mm_ebal()`: entropy balancing (Stata 11 required)
  * `mm_colvar()`: variance, by column
  * `mm_meancolvar()`: mean and variance, by column
  * `mm_variance0()`: population variance
  * `mm_meanvariance0()`: mean and population variance
  * `mm_mse()`: mean squared error
  * `mm_colmse()`: mean squared error, by column
  * `mm_sse()`: sum of squared errors
  * `mm_colsse()`: sum of squared errors, by column
  * `mm_benford()`: Benford distribution
  * `mm_cauchy()`: cumulative Cauchy-Lorentz dist.
  * `mm_cauchyden()`: Cauchy-Lorentz density
  * `mm_cauchytail()`: reverse cumulative Cauchy-Lorentz
  * `mm_invcauchy()`: inverse cumulative Cauchy-Lorentz
  * `mm_rbinomial()`: generate binomial random numbers
  * `mm_cebinomial()`: cond. expect. of binomial r.v.
  * `mm_root()`: Brent's univariate zero finder
  * `mm_nrroot()`: Newton-Raphson zero finder
  * `mm_minim()`: Brent's univariate minimum finder
  * `mm_finvert()`: univariate function inverter
  * `mm_integrate_sr()`: univariate function integration (Simpson's rule)
  * `mm_integrate_38()`: univariate function integration (Simpson's 3/8 rule)
  * `mm_ipolate()`: linear interpolation
  * `_mm_ipolate()`: linear interpolation (assuming sorted data)
  * `mm_fastipolate()`: linear interpolation (assuming sorted and unique data)
  * `mm_polint()`: polynomial inter-/extrapolation
  * `mm_sqrt()`: square root of a symmetric positive definite matrix
  * `mm_plot()`: Draw twoway plot
  * `_mm_plot()`: Draw twoway plot
  * `mm_group()`: create group index
  * `_mm_group()`: create group index, without sorting
  * `mm_panels()`: identify nested panel structure
  * `_mm_panels()`: identify panel sizes
  * `mm_npanels()`: identify number of panels
  * `mm_nunique()`: count number of unique values in vector
  * `mm_unique()`: obtain unique values from vector
  * `mm_unique_tag()`: tag unique values in vector
  * `mm_nuniqrows()`: count number of unique rows in matrix
  * `mm_uniqrows()`: obtain unique rows from matrix
  * `mm_uniqrows_tag()`: tag unique rows in matrix
  * `mm_diff()`: compute lagged differences
  * `mm_rowdiff()`: compute lagged differences within rows
  * `mm_coldiff()`: compute lagged differences within columns
  * `mm_isconstant()`: whether matrix is constant
  * `mm_issorted()`: whether vector is sorted
  * `mm_colrunsum()`: running sum of each column
  * `mm_prod()`: compute product of elements in matrix
  * `mm_rowprod()`: compute product within rows
  * `mm_colprod()`: compute product within columns
  * `mm_linbin()`: linear binning
  * `mm_fastlinbin()`: fast linear binning
  * `mm_exactbin()`: exact binning
  * `mm_makegrid()`: equally spaced grid points
  * `mm_seq()`: generate regular sequence
  * `mm_cut()`: categorize data vector
  * `mm_posof()`: find element in vector
  * `mm_which()`: positions of nonzero elements
  * `mm_locate()`: search an ordered vector
  * `mm_hunt()`: consecutive search
  * `mm_clip()`: clip/limit the values in a matrix
  * `mm_clipmin()`: limit the minimum
  * `mm_clipmax()`: limit the maximum
  * `mm_cond()`: matrix conditional operator
  * `mm_expand()`: duplicate single rows/columns
  * `_mm_expand()`: duplicate rows/columns in place
  * `mm_repeat()`: duplicate contents as a whole
  * `_mm_repeat()`: duplicate contents in place
  * `mm_sort()`: stable sorting
  * `mm_order()`: stable ordering
  * `mm_unorder2()`: stable version of `unorder()`
  * `mm_jumble2()`: stable version of `jumble()`
  * `mm__jumble2()`: stable version of `_jumble()`
  * `mm_pieces()`: break string into pieces
  * `mm_npieces()`: count number of pieces
  * `_mm_npieces()`: count number of pieces
  * `mm_regexr()`: regular expression replace
  * `mm_invtokens()`: reverse of `tokens()`
  * `mm_realofstr()`: convert string into real
  * `mm_strexpand()`: expand string argument
  * `mm_matlist()`: display a (real) matrix
  * `mm_insheet()`: read spreadsheet file
  * `mm_infile()`: read free-format file
  * `mm_outsheet()`: write spreadsheet file
  * `mm_callf()`: pass optional args to function
  * `mm_callf_setup()`: setup for `mm_callf()`

To install `moremata` from the SSC Archive, type

    . ssc install moremata, replace

in Stata. Stata version 9.2 or newer is required. Some functions may require
newer Stata versions.

---

Installation from GitHub:

    . net install moremata, replace from(https://raw.githubusercontent.com/benjann/moremata/master/)

---

Main changes:


    10aug2020
    - new mm_ddens() function for diffusion density estimation
    - mm_density():
      o ISJ bandwidth selector wrongly used grid size instead of number of obs when
        rescaling the bandwidth; this is fixed
      o increased padding of approximation grid to +/- 10% of data range (instead of
        +/- 5% percent)
      o D.n() now has an additional argument to set the padding proportion

    07aug2020
    - mm_density():
      o now using DPI if SJPI/ISJ fails
      o ISJ now uses same root-finding algorithm as SJPI
      o SJPI and DPI now compute the scale from the binned data
      o SJPI now uses min of sd and iqr as scale measure when computing the
        oversmoothed bandwidth; this is at odds with h_os(), but may add some
        robustness; furthermore, root finder now uses full precision
      o extension of automatic grid is now limited to 5% of range on either side
      o D.kernel() always selected gaussian; this is fixed
      o D.support(.,"",1) returned error; this is fixed

    06aug2020
    - new mm_density() funtion for (univariate) kernel density estimation
    - new mm_minim() funtion for univariate minimization without derivatives

    04aug2020
    - new mm_prod()/mm_rowprod()/mm_colprod() funtions to compute products of 
      elements in a matrix
    - new mm_seq() function to generate regular sequences

    17jul2020
    - improved quantile functions; underscore functions no longer assume weights
      to be nonzero and now allow multiple columns in P
    - new mm_issorted() function
    
    14jul2020
    - mm_quantile() has been rewritten; it now supports all 9 quantile definitions
      from from Hyndman and Fan (1996); weights are supported for all 
      definitions; new argument -fw- requests treating the weights as frequency
      weights
      ***
      IMPORTANT CHANGE:
          argument -altdef- in mm_quantile() and mm_iqrange() has been replaced
          by argument -def- that can take on values 1 to 9; altdef!=0 in the 
          previous version is equivalent to def=6 in the new version
      ***
    - new functions _mm_quantile(), _mm_median(), and _mm_iqrange() that assume
      sorted data
    - new functions mm_unique(), mm_unique_tag(), mm_uniqrows(), 
      mm_uniqrows_tag() to obtain or tag unique values in a vector or unique
      rows in a matrix; mm_uniqrows() differs from official uniqrows() in that
      it has an option to determin the order in which the result is returned
    - new functions _mm_nunique(), _mm_unique(), _mm_unique_tag(), 
      _mm_nuniqrows(), _mm_uniqrows(), and _mm_uniqrows_tag() to count, obtain,
      or tag unique values/rows without sorting the data
    - function mm_ipolate() is now faster, especially if there are ties
    - new _mm_ipolate() function that assumes sorted data
    - new mm_fastipolate() function that assumes sorted and unique data 
    - new mm_group() function for creating a group index
    - new mm_sort()/mm_order() functions for stable sorting
    - new mm_diff() function for lagged differences
    - new mm_clip() function to clip/limit values in a matrix
    - new mm_kderiv() function for kernel derivatives
    - new mm_ecdf2()/_mm_ecdf2() functions that return the CDF at unique values
      of X
    - argument -mid- in mm_ranks() did not make sense with ties=0 or ties=4; 
      this is fixed
    - function mm_relrank() has been reqritten; it now has additional 
      arguments support breaking ties and to compute nonnormalized ranks
    - new _mm_ecdf() function that assumes sorted data
    - new _mm_ranks() function that assumes sorted data
    - new _mm_relrank() function that assumes sorted data
    - mm_ranks() now uses quad precision in Stata 10 or newer
    - mm_ecdf(), mm_ranks(), and mm_relrank() now have separate help files
    - mm_colrunsum() now has argument -missing- to treat missing values as missing
      (instead of zero) and argument -quad- to request quad precision in Sata 10 
      or newer
    - mm_isconstant() now uses allof() instead of all() and is thus faster

    17apr2020
    - installation files added to GitHub distribution

    21aug2019
    - mm_ebal(): handling of collinearity/redundant constraints improved

    04may2019
    - strange problem caused by mm_ebal(): it left junk behind in memory; this
      had something to do with keeping an optimization object within a structure, 
      but then passing the structure as an argument to the optimization object; this
      fixed
    - mm_greedy() added

    26apr2019
    - mm_ebal() added

    30may2017
    - mm_sqrt() added

    01feb2017
    - mm_regexr() added

    01jun2015
    - mm_pieces() now supports unicode (Stata 14) 

    16may2014
    - mm_finvert() now has optional argument to pass on to &f()

    29jan2014
    - mm_integrate_sr() and mm_integrate_sr38() added

    19feb2009:
    - mm_collapse() added

    10feb2009:
    - mm_rbinomial(): note added that Stata 10.1 provides -rbinomial()-
    - mm_invtokens(): note added that Stata 10 provides -invtokens()-
    - new mm_pieces() functiom using genuine Mata code instead of extended macro 
      funtion -: piece-

    26mar2008:
    - mm_gini() updated so that it correctly handles ties. (Results depended on 
      sort order in case of ties)

    29feb2008:
    - mm_cond() added

    11jan2008
    - redirection of colrunsum in Stata 10 improved; _mm_colrunsum10() now faster
      if only 1 column
    - mm_invtokens() now also works with column vectors and has a -noclean- option
    - the default algorithm in mm_quantile() had precision problems if
      noninteger weights were specified
    - mm_quantile() now properly handles zero weights
    - mm_mgof() now displays progress dots
    - mm_mgof():
      - error message in cases where noninteger f is not allowed
      - mc method now rounds sum(f) to the nearest integer to prevent
        sampling (n-1) obs in case of imprecision

    29aug2007
    - mm_matlist() added
    - mm_colrunsum() now redirects itself to runningsum() if used in Stata 10

    07aug2007
    - mm_cauchy() functions added
    - mm_colrunsum() now no longer uses the mean update formula; the mean update formula
      is problematic with integers
    - mm_ranks() now has a normalize option (so that max(ecdf/relrank) is exactly 1)
      mm_gini(), mm_ecdf(), mm_relrank() updated
    - linbin/fastlinbin/exactbin now support data outside of grid
    - mm_ranks() has new syntax: new -mid- option for half-step (midpoint) method
      (replaces method==5); mm_relrank() now also has the mid option

    27jun2007
    - mm_benford() added
    - mm_upswor() now has a -nowarn- option
    - mm_ranks() changed (method=5 introduced; adjust removed; __mm_ranks() 
      for sorted data)
    - mm_relrank() now based on mm_ranks()
    - mm_nunique did not work with 'string rowvector' (because of transposeonly()); 
      this is fixed
    - mm_freq2() and _mm_freq2() added; _mm_freq() added
    - mm_freq() now allows matrix as input
    - mm_nuniqrows() is computed slightly differently now (faster if x has many 
      columns)
    - _mm_panels() is faster now
    - mm_isconstant() added
    - _mm_strexpand() added
    - bug with single quotes in strings with mm_pieces() fixed
    - mm_subset() etc. added
    - mm_mgof() added
    - mm_colrunsum(x) now works again if rows(x)==0 (the bug has been introduced
      on 12apr2007)
    - mm_which() now works if nothing is selected from a scalar

    12apr2007
    - mm_colrunsum() now uses the mean-update formula

    05apr2007
    - mm_pieces(), mm_npieces(), and _mm_pieces() added
    - mm_kern.mata: makes use of new capability of findexternal() to find
      functions; default kernel now epan2

    03aug2006
    - plot() added

    13jul2006
    - polint() added
    - kernel integrals for xK(x) and x^2K(x) added
    - slight changes to ipolate()
    - mm_nrroot, mm_finvert added
    - mm_locate, mm_hunt added
    - mm_root() added
    - fixed bug in mm_quantile
    - mm_kern: kernel functions added
    - default for m in makegrid() now 512 (previous: 401)
    - fixed bug with missings in variance0, mse, sse
    - w optional in quantile, iqrange, median, ecdf, ranks, freq, gini,
      histogram
    - P optional in quantile
    - g optional in histogram
    - mm_bs() and mm_jk() added
    - callf() added
    - slight change to mm_panels: info1 will be filled even if X1 is
      absent; info1 will contain two columns if Y==. or void
    - sse(), colsse() added
    - mse(), colmse() added
    - expand(), repeat() added
    - fw option deleted from linbin(), fastlinbin(), and exactbin()
    - nobs() added
    - nobs() now used in histogram()
    - quantile():
      * weighted version for altdef (only frequency weights)
      * speed improvements for unweighted algorithms
    - quantile() now has an altdef option (interpolation)
    - rank() now hat ties==4 option (order ties by w)
    - p in quantile(x,w,p) may now be matrix
    - q in relrank(x,w,q) may now be matrix

    23may2006
    - quantile, median, iqrange, ecdf, relrank, ranks now work with
      matrix X (statistics are computed for each column of X)
    - gini now works with matrix X (gini of each column of X)
    - fixed bug with with mm_sample() if stratified and n==0
    - mm_gini() added
    - mm_(mean)variance0(), mm_(mean)colvar() added
    - mm_rank() now has adjust option
    - mm_panels() etc: input now transmorphic vector
    - mm_nunique, mm_nuniqrows added
    - mm_ranks() added, mm_ecdf() now in terms of mm_ranks()
    - mm_npanels() added, mm_panels() can now be used with void strata
       and void cluster
    - mm_ipolate has new syntax (and is faster in most applications)
      (extrapolation not supported anymore; now using closest extremes)
    - mm_fastlinbin() added

    14apr2006
    - bug fixed in mm_sample() (nn[i] rather than n)
    - bug fixed in mm_sample() (strata[i,2] rather than cluster[i])
    - declarations fixed in rbinomial, cebinomial, outsheet
    - stable sort order in -mm_linbin()- and -mm_exactbin()-

    01apr2006
    - mm_unorder2(), mm_jumble2(), mm__jumble2() added
    - mm_sample() (etc.) added
    - mm_outsheet(): append/replace option
    - mm_panels() added
    - relrank(), ecdf(): range(1,I,1) changed to (1::I)
    - freq() added
    - cut() added
    - rbinomial() and cebinomial() added
    - posof() function added
    - insheet and infile:
      * much faster now (code based on cat() version 2)
      * now support reading specific range of file (line1-line2)

    15sep2005
    - released on SSC
