# moremata
Stata module providing various Mata functions

  * `mm_kern()`: various kernel functions
  * `mm_kint()`: kernel integral functions
  * `mm_kdel0()`: canonical bandwidth of kernel
  * `mm_quantile()`: quantile function
  * `mm_median()`: median
  * `mm_iqrange()`: inter-quartile range
  * `mm_ecdf()`: cumulative distribution function
  * `mm_relrank()`: grade transformation
  * `mm_ranks()`: ranks/cumulative frequencies
  * `mm_freq()`: compute frequency counts
  * `mm_histogram()`: produce histogram data
  * `mm_mgof()`: multinomial goodness-of-fit tests
  * `mm_collapse()`: summary statistics by subgroups
  * `_mm_collapse()`: summary statistics by subgroups
  * `mm_gini()`: Gini coefficient
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
  * `mm_finvert()`: univariate function inverter
  * `mm_integrate_sr()`: univariate function integration (Simpson's rule)
  * `mm_integrate_38()`: univariate function integration (Simpson's 3/8 rule)
  * `mm_ipolate()`: linear interpolation
  * `mm_polint()`: polynomial inter-/extrapolation
  * `mm_sqrt()`: square root of a symmetric positive definite matrix
  * `mm_plot()`: Draw twoway plot
  * `_mm_plot()`: Draw twoway plot
  * `mm_panels()`: identify nested panel structure
  * `_mm_panels()`: identify panel sizes
  * `mm_npanels()`: identify number of panels
  * `mm_nunique()`: count number of distinct values
  * `mm_nuniqrows()`: count number of unique rows
  * `mm_isconstant()`: whether matrix is constant
  * `mm_nobs()`: number of observations
  * `mm_colrunsum()`: running sum of each column
  * `mm_linbin()`: linear binning
  * `mm_fastlinbin()`: fast linear binning
  * `mm_exactbin()`: exact binning
  * `mm_makegrid()`: equally spaced grid points
  * `mm_cut()`: categorize data vector
  * `mm_posof()`: find element in vector
  * `mm_which()`: positions of nonzero elements
  * `mm_locate()`: search an ordered vector
  * `mm_hunt()`: consecutive search
  * `mm_cond()`: matrix conditional operator
  * `mm_expand()`: duplicate single rows/columns
  * `_mm_expand()`: duplicate rows/columns in place
  * `mm_repeat()`: duplicate contents as a whole
  * `_mm_repeat()`: duplicate contents in place
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

To install moremata in Stata, type

    . ssc install moremata, replace

or download `moremata.zip` from
[RePEc](http://ideas.repec.org/c/boc/bocode/s455001.html)
and follow the manual installation instructions in the readme therein.

Stata version 9.2 or newer is required.
