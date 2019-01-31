*! version 1.0.0, Ben Jann, 19jun2007
version 9.2
mata:

real vector mm_benford(
 real vector digit0,
 | real scalar pos0,
   real scalar base0)
{
    real scalar    pos, base
    real vector    digit

// check arguments and set defaults
    if (args()<3 | base0>=.) base = 10
    else {
        base = trunc(base0)
        if (base <0 | base >10) _error(3300, "base must be in [2,10]")
    }
    if (args()<2 | pos0>=.) pos = 1
    else {
        pos = trunc(pos0)
        if (pos<1) _error(3300, "pos must be 1 or larger")
    }
    digit = trunc(digit0)
    if (any(digit:<0) | any(digit:>=base))
      _error(3300, "digit must be in [0,base-1] (base is "+strofreal(base)+")")

// case 1: pos==1 (simple formula)
    if (pos==1) return(ln(1 :+ 1:/digit)/ln(base))

// case 2: pos > 1
    if (cols(digit)==1) return( rowsum(ln(1 :+ 1:/(
     J(rows(digit),1,base)*(base^(pos-2)..base^(pos-1)-1) :+ digit))) / ln(base))
    else                return( colsum(ln(1 :+ 1:/(
     (base^(pos-2)::base^(pos-1)-1)*J(1,cols(digit),base) :+ digit))) / ln(base))
}

end
