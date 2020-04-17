*! version 1.0.1, Ben Jann, 16jun2006
version 9.2

local o1 "o1"
local po1 "*p.o1"
forv i=2/10 {
	local o`i' "`o`=`i'-1'', o`i'"
	local po`i' "`po`=`i'-1'', *p.o`i'"
}

mata:

struct mm_callf_o10 {
	pointer(function) scalar f
	real scalar              n
	pointer scalar           `o10'
}

struct mm_callf_o10 scalar mm_callf_setup(
 pointer(function) scalar f, real scalar n, | `o10')
{
	struct mm_callf_o10 scalar p

	p.f  = f
	p.n  = n
	p.o1 = &o1 ; p.o2 = &o2 ; p.o3 = &o3 ; p.o4 = &o4 ; p.o5 = &o5
	p.o6 = &o6 ; p.o7 = &o7 ; p.o8 = &o8 ; p.o9 = &o9 ; p.o10 = &o10
	return(p)
}

transmorphic mm_callf(struct mm_callf_o10 p, | `o5')
{
	if (args()==1) return(_mm_callf0(p))
	if (args()==2) return(_mm_callf1(p, `o1'))
	if (args()==3) return(_mm_callf2(p, `o2'))
	if (args()==4) return(_mm_callf3(p, `o3'))
	if (args()==5) return(_mm_callf4(p, `o4'))
	               return(_mm_callf5(p, `o5'))
}

transmorphic _mm_callf0(struct mm_callf_o10 p)
{
	if (p.n<=0) return( (*p.f)() )
	if (p.n==1) return( (*p.f)(`po1') )
	if (p.n==2) return( (*p.f)(`po2') )
	if (p.n==3) return( (*p.f)(`po3') )
	if (p.n==4) return( (*p.f)(`po4') )
	if (p.n==5) return( (*p.f)(`po5') )
	if (p.n==6) return( (*p.f)(`po6') )
	if (p.n==7) return( (*p.f)(`po7') )
	if (p.n==8) return( (*p.f)(`po8') )
	if (p.n==9) return( (*p.f)(`po9') )
	            return( (*p.f)(`po10') )
}

transmorphic _mm_callf1(struct mm_callf_o10 p, `o1')
{
	if (p.n<=0) return( (*p.f)(`o1') )
	if (p.n==1) return( (*p.f)(`o1', `po1') )
	if (p.n==2) return( (*p.f)(`o1', `po2') )
	if (p.n==3) return( (*p.f)(`o1', `po3') )
	if (p.n==4) return( (*p.f)(`o1', `po4') )
	if (p.n==5) return( (*p.f)(`o1', `po5') )
	if (p.n==6) return( (*p.f)(`o1', `po6') )
	if (p.n==7) return( (*p.f)(`o1', `po7') )
	if (p.n==8) return( (*p.f)(`o1', `po8') )
	if (p.n==9) return( (*p.f)(`o1', `po9') )
	            return( (*p.f)(`o1', `po10') )
}

transmorphic _mm_callf2(struct mm_callf_o10 p, `o2')
{
	if (p.n<=0) return( (*p.f)(`o2') )
	if (p.n==1) return( (*p.f)(`o2', `po1') )
	if (p.n==2) return( (*p.f)(`o2', `po2') )
	if (p.n==3) return( (*p.f)(`o2', `po3') )
	if (p.n==4) return( (*p.f)(`o2', `po4') )
	if (p.n==5) return( (*p.f)(`o2', `po5') )
	if (p.n==6) return( (*p.f)(`o2', `po6') )
	if (p.n==7) return( (*p.f)(`o2', `po7') )
	if (p.n==8) return( (*p.f)(`o2', `po8') )
	if (p.n==9) return( (*p.f)(`o2', `po9') )
	            return( (*p.f)(`o2', `po10') )
}

transmorphic _mm_callf3(struct mm_callf_o10 p, `o3')
{
	if (p.n<=0) return( (*p.f)(`o3') )
	if (p.n==1) return( (*p.f)(`o3', `po1') )
	if (p.n==2) return( (*p.f)(`o3', `po2') )
	if (p.n==3) return( (*p.f)(`o3', `po3') )
	if (p.n==4) return( (*p.f)(`o3', `po4') )
	if (p.n==5) return( (*p.f)(`o3', `po5') )
	if (p.n==6) return( (*p.f)(`o3', `po6') )
	if (p.n==7) return( (*p.f)(`o3', `po7') )
	if (p.n==8) return( (*p.f)(`o3', `po8') )
	if (p.n==9) return( (*p.f)(`o3', `po9') )
	            return( (*p.f)(`o3', `po10') )
}

transmorphic _mm_callf4(struct mm_callf_o10 p, `o4')
{
	if (p.n<=0) return( (*p.f)(`o4') )
	if (p.n==1) return( (*p.f)(`o4', `po1') )
	if (p.n==2) return( (*p.f)(`o4', `po2') )
	if (p.n==3) return( (*p.f)(`o4', `po3') )
	if (p.n==4) return( (*p.f)(`o4', `po4') )
	if (p.n==5) return( (*p.f)(`o4', `po5') )
	if (p.n==6) return( (*p.f)(`o4', `po6') )
	if (p.n==7) return( (*p.f)(`o4', `po7') )
	if (p.n==8) return( (*p.f)(`o4', `po8') )
	if (p.n==9) return( (*p.f)(`o4', `po9') )
	            return( (*p.f)(`o4', `po10') )
}

transmorphic _mm_callf5(struct mm_callf_o10 p, `o5')
{
	if (p.n<=0) return( (*p.f)(`o5') )
	if (p.n==1) return( (*p.f)(`o5', `po1') )
	if (p.n==2) return( (*p.f)(`o5', `po2') )
	if (p.n==3) return( (*p.f)(`o5', `po3') )
	if (p.n==4) return( (*p.f)(`o5', `po4') )
	if (p.n==5) return( (*p.f)(`o5', `po5') )
	if (p.n==6) return( (*p.f)(`o5', `po6') )
	if (p.n==7) return( (*p.f)(`o5', `po7') )
	if (p.n==8) return( (*p.f)(`o5', `po8') )
	if (p.n==9) return( (*p.f)(`o5', `po9') )
	            return( (*p.f)(`o5', `po10') )
}

end
