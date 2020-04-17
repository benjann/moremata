*! version 1.0.1, Ben Jann, 30aug2006
version 9.2

local kernels ///
 epanechnikov ///
 epan2 ///
 biweight ///
 triweight ///
 cosine ///
 gaussian ///
 parzen ///
 rectangle ///
 triangle

local klist tokens("`kernels'")

local kname `"(k!="" ? k : "epan2")"'

mata:

string scalar _mm_unabkern(|string scalar k)
{
	return(mm_strexpand(k, `klist', "epan2", 0, k+": invalid kernel"))
}


pointer scalar _mm_findkern(|string scalar k)
{
	pointer(real matrix function) scalar K
	K = findexternal("mm_kern_"+`kname'+"()")
	if (K==NULL) _error(3499,"mm_kern_"+`kname'+"() not found")
	return(K)
}

pointer scalar _mm_findkint(|string scalar k)
{
	pointer(real matrix function) scalar K
	K = findexternal("mm_kint_"+`kname'+"()")
	if (K==NULL) _error(3499,"mm_kint_"+`kname'+"() not found")
	return(K)
}

pointer scalar _mm_findkdel0(|string scalar k)
{
	pointer(real scalar function) scalar K
	K = findexternal("mm_kdel0_"+`kname'+"()")
	if (K==NULL) _error(3499,"mm_kdel0_"+`kname'+"() not found")
	return(K)
}

real matrix mm_kern(string scalar k, real matrix x)
{
	return((*_mm_findkern(_mm_unabkern(k)))(x))
}

real matrix mm_kint(string scalar k, real scalar r, |real matrix x)
{
	return(mm_callf(mm_callf_setup(_mm_findkint(_mm_unabkern(k)), args()-2, x), r))
}

real matrix mm_kdel0(string scalar k)
{
	return((*_mm_findkdel0(_mm_unabkern(k)))())
}

// mm_kern_*: kernel function
// mm_kint_*(1,x): integral of K(z) from -infty to x
// mm_kint_*(2,x): integral of K(z)^2 from -infty to x
// mm_kint_*(2): kernel "roughness"
// mm_kint_*(3,x): integral of z*K(z) from -infty to x
// mm_kint_*(4,x): integral of z^2*K(z) from -infty to x
// mm_kint_*(4): kernel variance
// mm_kvar_*: variance
// mm_kdel0_*: canonical bandwidth

// (scaled) Epanechnikov kernel function
real matrix mm_kern_epanechnikov(real matrix x)
{
	return((.75:-.15*x:^2)/sqrt(5):*(abs(x):<sqrt(5)))
}
real matrix mm_kint_epanechnikov(real scalar r, |real matrix x)
{
	if (r==1 & args()<2) return(1)
	if (r==1)
	 return((.5:+(.75*x-.05*x:^3)/sqrt(5)):*(abs(x):<sqrt(5)) + (x:>=sqrt(5)))
	if (r==2 & args()<2) return(3/5^1.5)
	if (r==2)
	 return((.3/sqrt(5):+(.1125*x-.015*x:^3+.0009*x:^5)):*(abs(x):<sqrt(5)) + 3/5^1.5*(x:>=sqrt(5)))
	if (r==3 & args()<2) return(0)
	if (r==3)
	 return((-.1875*sqrt(5):+.375/sqrt(5)*x:^2-.0375/sqrt(5)*x:^4):*(abs(x):<sqrt(5)))
	if (r==4 & args()<2) return(1)
	if (r==4)
	 return((.5:+.25/sqrt(5)*x:^3-.03/sqrt(5)*x:^5):*(abs(x):<sqrt(5)) + (x:>=sqrt(5)))
	else _error(3300)
}
//real scalar mm_kvar_epanechnikov() return(1)
real scalar mm_kdel0_epanechnikov() return((3/5^1.5)^.2)


// Epanechnikov kernel function
real matrix mm_kern_epan2(real matrix x)
{
	return((.75:-.75*x:^2):*(abs(x):<1))
}
real matrix mm_kint_epan2(real scalar r, |real matrix x)
{
	if (r==1 & args()<2) return(1)
	if (r==1)
	 return((.5:+.75*x-.25*x:^3):*(abs(x):<1) + (x:>=1))
	if (r==2 & args()<2) return(.6)
	if (r==2)
	 return((.3:+.5625*x-.375*x:^3+.1125*x:^5):*(abs(x):<1) + .6*(x:>=1))
	if (r==3 & args()<2) return(0)
	if (r==3)
	 return((-.1875:+.375*x:^2-.1875*x:^4):*(abs(x):<1))
	if (r==4 & args()<2) return(.2) // 1/5
	if (r==4)
	 return((.1:+.25*x:^3-.15*x:^5):*(abs(x):<1) + .2*(x:>=1))
	else _error(3300)
}
//real scalar mm_kvar_epan2() return(.2)
real scalar mm_kdel0_epan2() return(15^.2)


// Biweight kernel function
real matrix mm_kern_biweight(real matrix x)
{
	return(.9375*(1:-x:^2):^2:*(abs(x):<1))
}
real matrix mm_kint_biweight(| real scalar r, real matrix x)
{
	if (r==1 & args()<2) return(1)
	if (r==1)
	 return((.5:+.9375*(x-2/3*x:^3+.2*x:^5)):*(abs(x):<1) + (x:>=1))
	if (r==2 & args()<2) return(5/7)
	if (r==2)
	 return((5/14:+225/256*x-75/64*x:^3+135/128*x:^5-225/448*x:^7+25/256*x:^9):*(abs(x):<1) + 5/7*(x:>=1))
	if (r==3 & args()<2) return(0)
	if (r==3)
	 return((-.15625:+.46875*x:^2-.46875*x:^4 +.15625*x:^6):*(abs(x):<1))
	if (r==4 & args()<2) return(1/7)
	if (r==4)
	 return((1/14:+.3125*x:^3-.375*x:^5+15/112*x:^7):*(abs(x):<1) + (x:>=1)/7)
	else _error(3300)
}
//real scalar mm_kvar_biweight() return(1/7)
real scalar mm_kdel0_biweight() return(35^.2)


// Triweight kernel function
real matrix mm_kern_triweight(real matrix x)
{
	return(1.09375*(1:-x:^2):^3:*(abs(x):<1))
}
real matrix mm_kint_triweight(| real scalar r, real matrix x)
{
	if (r==1 & args()<2) return(1)
	if (r==1)
	 return((.5:+1.09375*(x-x:^3+.6*x:^5-1/7*x:^7)):*(abs(x):<1) + (x:>=1))
	if (r==2 & args()<2) return(350/429)
	if (r==2)
	 return((175/429:+1225/1024*(x-2*x:^3+3*x:^5-20/7*x:^7+5/3*x:^9-6/11*x:^11+1/13*x:^13)):*(abs(x):<1) + 350/429*(x:>=1))
	if (r==3 & args()<2) return(0)
	if (r==3)
	 return((-35/256:+35/64*x:^2-105/128*x:^4+35/64*x:^6-35/256*x:^8):*(abs(x):<1))
	if (r==4 & args()<2) return(1/9)
	if (r==4)
	 return((1/18:+35/96*x:^3-21/32*x:^5+15/32*x:^7-35/288*x:^9):*(abs(x):<1) + (x:>=1)/9)
	else _error(3300)
}
//real scalar mm_kvar_triweight() return(1/9)
real scalar mm_kdel0_triweight() return((9450/143)^.2)


// Cosine trace
real matrix mm_kern_cosine(real matrix x)
{
	return((1:+cos(2*pi()*x)):*(abs(x):<.5))
}
real matrix mm_kint_cosine(| real scalar r, real matrix x)
{
	if (r==1 & args()<2) return(1)
	if (r==1)
	 return((.5:+x+sin(2*pi()*x)/(2*pi())):*(abs(x):<.5) + (x:>=.5))
	if (r==2 & args()<2) return(1.5)
	if (r==2)
	 return((.75:+1.5*x+sin(2*pi()*x)/pi()+cos(2*pi()*x):*sin(2*pi()*x)/(4*pi())):*(abs(x):<.5) + 1.5*(x:>=.5))
	if (r==3 & args()<2) return(0)
	if (r==3)
	 return((-.125+.25/pi()^2:+cos(2*pi()*x)/(4*pi()^2)+x:*sin(2*pi()*x)/(2*pi())+.5*x:^2):*(abs(x):<.5))
	if (r==4 & args()<2) return(.5*(1/6-1/pi()^2))
	if (r==4)
	 return((1/24-.25/pi()^2:-sin(2*pi()*x)/(4*pi()^3)+x:*cos(2*pi()*x)/(2*pi()^2)+x:^2:*sin(2*pi()*x)/(2*pi())+x:^3/3):*(abs(x):<.5) + .5*(1/6-1/pi()^2)*(x:>=.5))
	else _error(3300)
}
//real scalar mm_kvar_cosine() return(.5*(1/6-1/pi()^2))
real scalar mm_kdel0_cosine() return((3/(.5*(1/6-1/pi()^2)^2))^.2)


// Gaussian kernel function
real matrix mm_kern_gaussian(real matrix x)
{
	return(normalden(x)) // alternatively: exp(-0.5*((x)^2))/sqrt(2*pi())
}
real matrix mm_kint_gaussian(| real scalar r, real matrix x)
{
	if (r==1 & args()<2) return(1)
	if (r==1)            return(normal(x))
	if (r==2 & args()<2) return(.5/sqrt(pi())) // = normalden(0,sqrt(2))
	if (r==2)            return(.5/sqrt(pi()) * normal(sqrt(2)*x))
	if (r==3 & args()<2) return(0)
	if (r==3)            return(-1/sqrt(2*pi())*exp(-.5*x:^2))
	if (r==4 & args()<2) return(1)
	if (r==4)            return(-1/sqrt(2*pi())*x:*exp(-.5*x:^2) + normal(x))
	else _error(3300)
}
//real scalar mm_kvar_gaussian() return(1)
real scalar mm_kdel0_gaussian() return((1/(4*pi()))^.1)


// Parzen kernel function
real matrix mm_kern_parzen(real matrix x)
{
	return(((4/3:-8*x:^2+8*abs(x):^3):*(abs(x):<=.5)) +
	 ((8*(1:-abs(x)):^3/3):*(abs(x):>1/2:&abs(x):<=1)))
}
real matrix mm_kint_parzen(| real scalar r, real matrix x)
{
	if (r==1 & args()<2) return(1)
	if (r==1)
	 return((2/3:+8/3*(x+1.5*x:^2+x:^3+.25*x:^4)):*(x:>=-1:&x:<-.5) +
	  (.5:+4/3*x-8/3*x:^3-2*x:^4):*(x:>=-.5:&x:<0) +
	  (.5:+4/3*x-8/3*x:^3+2*x:^4):*(x:>=0:&x:<=.5) +
	  (1/3:+8/3*(x-1.5*x:^2+x:^3-.25*x:^4)):*(x:>.5:&x:<=1) + (x:>1))
	if (r==2 & args()<2) return(302/315)
	if (r==2)
	 return((64/63:+64/9*(x+3*x:^2+5*x:^3+5*x:^4+3*x:^5+x:^6+x:^7/7)):*(x:>=-1:&x:<-.5) +
	  (151/315:+16/9*(x-4*x:^3-3*x:^4+36/5*x:^5+12*x:^6+36/7*x:^7)):*(x:>=-.5:&x:<0) +
	  (151/315:+16/9*(x-4*x:^3+3*x:^4+36/5*x:^5-12*x:^6+36/7*x:^7)):*(x:>=0:&x:<=.5) +
	  (-2/35:+64/9*(x-3*x:^2+5*x:^3-5*x:^4+3*x:^5-x:^6+x:^7/7)):*(x:>.5:&x:<=1) + 302/315*(x:>1))
	if (r==3 & args()<2) return(0)
	if (r==3)
	 return((-2/15:+4/3*x:^2+8/3*x:^3+2*x:^4+8/15*x:^5):*(x:>=-1:&x:<-.5) +
	  (-7/60:+2/3*x:^2-2*x:^4-1.6*x:^5):*(x:>=-.5:&x:<0) +
	  (-7/60:+2/3*x:^2-2*x:^4+1.6*x:^5):*(x:>=0:&x:<=.5) +
	  (-2/15:+4/3*x:^2-8/3*x:^3+2*x:^4-8/15*x:^5):*(x:>.5:&x:<=1))
	if (r==4 & args()<2) return(1/12)
	if (r==4)
	 return((2/45:+8/9*x:^3+2*x:^4+1.6*x:^5+4/9*x:^6):*(x:>=-1:&x:<-.5) +
	  (1/24:+4/9*x:^3-1.6*x:^5-4/3*x:^6):*(x:>=-.5:&x:<0) +
	  (1/24:+4/9*x:^3-1.6*x:^5+4/3*x:^6):*(x:>=0:&x:<=.5) +
	  (7/180:+8/9*x:^3-2*x:^4+1.6*x:^5-4/9*x:^6):*(x:>.5:&x:<=1) + (x:>1)/12)
	else _error(3300)
}
//real scalar mm_kvar_parzen() return(1/12)
real scalar mm_kdel0_parzen() return((4832/35)^.2)

// Rectangle kernel function
real matrix mm_kern_rectangle(real matrix x)
{
	return(.5*(round(abs(x),1e-8):<1))
} // 1e-8 same as in kdensity.ado, version 2.3.6 26jun2000 (Stata 7)
real matrix mm_kint_rectangle(| real scalar r, real matrix x)
{
	if (r==1 & args()<2) return(1)
	if (r==1)
	 return((.5:+.5*x):*(round(abs(x),1e-8):<1) + (round(x,1e-8):>=1))
	if (r==2 & args()<2) return(.5)
	if (r==2)
	 return((.25:+.25*x):*(round(abs(x),1e-8):<1) + .5*(round(x,1e-8):>=1))
	if (r==3 & args()<2) return(0)
	if (r==3)
	 return((-.25:+.25*x:^2):*(round(abs(x),1e-8):<1))
	if (r==4 & args()<2) return(1/3)
	if (r==4)
	 return((1/6:+x:^3/6):*(round(abs(x),1e-8):<1) + (round(x,1e-8):>=1)/3)
	else _error(3300)
}
//real scalar mm_kvar_rectangle() return(1/3)
real scalar mm_kdel0_rectangle() return((9/2)^.2)


// Triangle kernel function
real matrix mm_kern_triangle(real matrix x)
{
	return((1:-abs(x)):*(abs(x):<1))
}
real matrix mm_kint_triangle(| real scalar r, real matrix x)
{
	if (r==1 & args()<2) return(1)
	if (r==1)
	 return((.5:+x+.5*x:^2):*(x:>-1:&x:<0) +
	  (.5:+x-.5*x:^2):*(x:>=0:&x:<1) + (x:>=1))
	if (r==2 & args()<2) return(2/3)
	if (r==2)
	 return((x+x:^2+(1:+x:^3)/3):*(x:>-1:&x:<0) +
	  (x-x:^2+(1:+x:^3)/3):*(x:>=0:&x:<1) + 2/3*(x:>=1))
	if (r==3 & args()<2) return(0)
	if (r==3)
	 return((-1/6:+.5*x:^2+x:^3/3):*(x:>-1:&x:<0) +
	  (-1/6:+.5*x:^2-x:^3/3):*(x:>=0:&x:<1))
	if (r==4 & args()<2) return(1/6)
	if (r==4)
	 return((1/12:+x:^3/3+.25*x:^4):*(x:>-1:&x:<0) +
	  (1/12:+x:^3/3-.25*x:^4):*(x:>=0:&x:<1) + (x:>=1)/6)
	else _error(3300)
}
//real scalar mm_kvar_triangle() return(1/6)
real scalar mm_kdel0_triangle() return(24^.2)

end
