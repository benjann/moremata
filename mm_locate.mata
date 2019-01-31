*! version 1.0.0, Ben Jann, 07jul2006

version 9.2

mata:

void mm_locate(
// translation of -locate- from Press et al. (1992), Numerical
// Recipes in C, p. 117, http://www.numerical-recipes.com
// Description: "Given an array xx[1..n], and given a value x,
// returns a value j such that x is between xx[j] and xx[j+1].
// xx must be monotonic, either increasing or decreasing. j=0
// or j=n is returned to indicate that x is out of range."
 numeric vector xx,
 numeric scalar x,
 real scalar j)
{
	real scalar  ju, jm, jl
	real scalar  ascnd
	real scalar  n

	n = length(xx)
	jl = 0
	ju = n + 1
	ascnd = (xx[n] >= xx[1])
	while (ju-jl > 1) {
		jm = trunc((ju+jl)/2)
		if (x >= xx[jm] == ascnd)
		  jl = jm
		else
		  ju = jm
	}
	if (x == xx[1])     j = 1
	else if(x == xx[n]) j = n - 1
	else                j = jl
}

void mm_hunt(
// translation of -hunt- from Press et al. (1992), Numerical
// Recipes in C, p. 118-9, http://www.numerical-recipes.com
// Description: "Given an array xx[1..n], and given a value x,
// returns a value jlo such that x is between xx[jlo] and
// xx[jlo+1]. xx[1..n] must be monotonic, either increasing or
// decreasing. jlo=0 or jlo=n is returned to indicate that x is
// out of range. jlo on input is taken as the initial guess for
// jlo on output.
 numeric vector xx,
 numeric scalar x,
 real scalar jlo)
{
	real scalar  jm, jhi, inc
	real scalar  ascnd
	real scalar  n

	n = length(xx)
	ascnd=(xx[n] >= xx[1])
	if (jlo <= 0 | jlo > n) {
		jlo = 0
		jhi = n + 1
	}
	else {
		inc = 1
		if (x >= xx[jlo] == ascnd) {
			if (jlo == n) return
			jhi = jlo + 1
			while (x >= xx[jhi] == ascnd) {
				jlo = jhi
				inc = inc + inc
				jhi = jlo + inc
				if (jhi > n) {
					jhi = n + 1
					break
				}
			}
		}
		else {
			if (jlo == 1) {
				jlo = 0
				return
			}
			jhi = jlo--
			while (x < xx[jlo] == ascnd) {
				jhi = jlo
				inc = inc + inc
				if (inc >= jhi) {
					jlo = 0
					break
				}
				else  jlo = jhi - inc
			}
		}
	}
	while (jhi-jlo != 1) {
		jm = trunc((jhi+jlo)/2)
		if (x >= xx[jm] == ascnd)
		 jlo = jm
		else
		 jhi = jm
	}
	if (x == xx[n]) jlo = n - 1
	if (x == xx[1]) jlo = 1
}

end
