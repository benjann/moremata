*! version 2.0.1, Ben Jann, 18mar2006
*! based on cat(), version 2.0.0  23jan2006
version 9.0
mata:

string matrix mm_insheet(string scalar filename,
 | string scalar del, real scalar uline1, real scalar line2)
{
	real scalar    fh, fpos, line1, i, r, c, j, delpos, ji
	string matrix  EOF
	string matrix  res
	string scalar  line

	/* ------------------------------------------------------------ */
	/* setup */
	if (del=="") del = char(9)
	fh  = fopen(filename, "r")
	EOF = J(0, 0, "")
	line1 = floor(uline1)
	if (line1<1 | line1>=.) line1 = 1

	/* ------------------------------------------------------------ */
	/* nothing to read case  */
	if (line1>line2) {
		fclose(fh)
		return(J(0, 0, ""))
	}

	/* ------------------------------------------------------------ */
	 /* skip opening lines  */
	for (i=1; i<line1; i++) {
		if (fget(fh)==EOF) {
			fclose(fh)
			return(J(0, 0, ""))
		}
	}
	fpos = ftell(fh)

	/* ------------------------------------------------------------ */
	/* count lines and columns  */
	j = 0
	for (i=line1; i<=line2; i++) {
		if ((line=fget(fh))==EOF) break
		ji = 1
		while (delpos = strpos(line, del)) {
			ji++
			line = substr(line, delpos+1, .)
		}
		if (ji>j) j = ji
	}
	res = J(r = i-line1, c = j, "")

	/* ------------------------------------------------------------ */
	/* read lines */
	fseek(fh, fpos, -1)
	for (i=1; i<=r; i++) {
		if ((line=fget(fh))==EOF) {
			/* unexpected EOF -- file must have changed */
			fclose(fh)
			if (--i) return(res[|1,1 \ i,c|])
			return(J(0,c,""))
		}
		for (j=1; j<=c; j++) {
			delpos = strpos(line, del)
			if (delpos==0) {
				res[i,j] = line
				break
			}
			res[i,j] = substr(line, 1, delpos-1)
			line = substr(line, delpos+1, .)
		}
	}
	fclose(fh)
	return(res)
}

end
