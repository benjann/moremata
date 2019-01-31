*! version 1.0.4  31aug2007  Ben Jann
version 9.2
mata:

void mm_matlist(
    real matrix X,
    | string matrix fmt,
      real scalar b,
      string vector rstripe,
      string vector cstripe,
      string scalar rtitle,
      transmorphic matrix res
    )
{
    real scalar             i, j, rw, lw, j0, j1, l, r, save
    real rowvector          wd
    string rowvector        sfmt
    pointer(real scalar) scalar     fi, fj
    pointer(string vector) scalar   rstr, cstr

    if (args()<2) fmt = "%9.0g"
    if (args()<3) b = 1
    save = (args()==7)
    //wd = strlen(sprintf(fmt,-1/3))

    if ((rows(fmt)!=1 & rows(fmt)!=rows(X)) |
        (cols(fmt)!=1 & cols(fmt)!=cols(X))) _error(3200)
    if (!anyof((0,1,2,3),b)) _error(3300)

    fi = (rows(fmt)!=1 ? &i : &1)
    fj = (cols(fmt)!=1 ? &j : &1)

    wd = J(1,cols(X),0)
    for (i=1;i<=rows(X);i++) {
        for (j=1;j<=cols(X);j++) {
            wd[j] = max((wd[j],strlen(sprintf(fmt[*fi,*fj],X[i,j]))))
        }
    }

    if (length(X)==0) return

    if (length(rstripe)<1)  rstr = &strofreal(1::rows(X))
    else rstr = &rstripe
    if (length(cstripe)<1)  cstr = &strofreal(1..cols(X))
    else if (rows(cstripe)!=1)  cstr = &(cstripe')
    else                        cstr = &cstripe

    if ((*rstr!="" & length(*rstr)!=rows(X)) |
        (*cstr!="" & length(*cstr)!=cols(X))) _error(3200)

    if (*rstr!="")  rw = max((max(strlen(*rstr)),strlen(rtitle)))
    else            rw = -1
    if (*cstr!="")  wd = colmax((strlen(*cstr)\wd)) :+ 2
    else            wd = wd :+ 2

    sfmt = "%" :+ strofreal(wd) :+ "s"

    if (save) {
        _mm_matlist(X, fmt, b, *rstr, *cstr, rtitle,
            wd, rw, sfmt, 1, cols(X), 1, 1, res)
        return
    }

    lw = st_numscalar("c(linesize)") - ((2+rw+2*(b>0))+2*(b==1))
    if ((sum(wd)+cols(X))<=lw) {
        _mm_matlist(X, fmt, b, *rstr, *cstr, rtitle,
            wd, rw, sfmt, 1, cols(X), 1, 1)
        return
    }

    j0 = 1
    l = 1
    r = 0
    while (j0<=cols(X)) {
        j1 = j0
        while (sum(wd[|j0\j1|]:+1)<=lw) {
            if (++j1>cols(X)) {
                r = 1
                break
            }
        }
        j1 = max((j1-1,j0))
        _mm_matlist(X, fmt, b, *rstr, *cstr, rtitle,
            wd, rw, sfmt, j0, j1, l, r)
        l = 0
        j0 = j1+1
    }
}

string colvector mm_smatlist(
    real matrix X,
    | string matrix fmt,
      real scalar b,
      string vector rstripe,
      string vector cstripe,
      string scalar rtitle
    )
{
    string colvector res

    if (args()<2) fmt = "%9.0g"
    if (args()<3) b = 1
    if (args()<4) rstripe = J(0,1,"")
    if (args()<5) cstripe = J(1,0,"")
    if (args()<5) rtitle = ""
    mm_matlist(X, fmt, b, rstripe, cstripe, rtitle, res)
    return(res)
}

void _mm_matlist(
    real matrix X,
    string matrix fmt,
    real scalar b,
    string vector rstripe,
    string vector cstripe,
    string scalar rtitle,
    real rowvector wd,
    real scalar rw,
    string rowvector sfmt,
    real scalar j0,
    real scalar j1,
    real scalar l,
    real scalar r,
    | transmorphic matrix res
    )
{
    real scalar          i, j, di, ri
    string scalar        rfmt, line
    pointer(real scalar) scalar fi, fj

    di = args()<14
    fi = (rows(fmt)!=1 ? &i : &1)
    fj = (cols(fmt)!=1 ? &j : &1)
    rfmt = "%"+strofreal(rw)+"s"

    if (di==0) {
        ri = 0
        res = J(rows(X)+(1+(b==3))*(cstripe!="")+(b>0)+(b==1)
            +(b==3)+(l==0&(b==0|b==2)),1,"")
    }

    if (l==0 & (b==0 | b==2)) {
        if (di) printf("\n")
        else res[++ri] = ""
    }
    if (cstripe!="") {
        if (b==3) {
            line = sprintf("{txt}" + "{hline " + strofreal(2+rw+1) +  "}{c TT}"
                + "{hline " + strofreal(sum(wd[|j0\j1|]:+1)-1) + "}")
            if (di) display(line)
            else res[++ri] = line
        }
        line = sprintf("{txt}  ")
        if (rstripe!="") line = line + sprintf(rfmt+" ", rtitle)
        if (b>0) line = line + sprintf((b==1 ? " " : "{c |}"))
        for (j=j0;j<=j1;j++) {
            line = line + sprintf(sfmt[j], cstripe[j])
            if (j<j1) line = line + sprintf(" ")
        }
        if (di) display(line)
        else res[++ri] = line
    }
    if (b) {
        line = sprintf("{txt}" + (b==1 ? (2+rw+1)*" " +
            (l ? "{c TLC}" : " ") : "{hline " + strofreal(2+rw+1) + "}{c +}")
            + "{hline " + strofreal(sum(wd[|j0\j1|]:+1)+1-2*(b!=1)) + "}" +
            (b==1 ? (r ? "{c TRC}" : "") : ""))
        if (di) display(line)
        else res[++ri] = line
    }
    for (i=1;i<=rows(X);i++) {
        line = sprintf("{txt}  ")
        if (rstripe!="") line = line + sprintf(rfmt+" ", rstripe[i])
        if (b>0) line = line + sprintf((l | b!=1 ? "{c |}" : " "))
        line = line + sprintf("{res}")
        for (j=j0;j<=j1;j++) {
            line = line + sprintf(sfmt[j],sprintf(fmt[*fi,*fj], X[i,j]))
            if (j<j1) line = line + sprintf(" ")
        }
        if (b==1 & r) line = line + sprintf("  {txt}{c |}")
        if (di) display(line)
        else res[++ri] = line
    }
    if (b==1|b==3) {
        line = sprintf("{txt}" + (b==1 ? (2+rw+1)*" " +
            (l ? "{c BLC}" : " ") : "{hline " + strofreal(2+rw+1) + "}{c BT}")
            + "{hline " + strofreal(sum(wd[|j0\j1|]:+1)+1-2*(b!=1)) + "}" +
            (b==1 ? (r ? "{c BRC}" : "") : ""))
        if (di) display(line)
        else res[++ri] = line
    }
}

end
