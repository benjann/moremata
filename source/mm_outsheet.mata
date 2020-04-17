*! version 1.0.4, Ben Jann, 14apr2006
version 9.0
mata:

void mm_outsheet(string scalar fn, string matrix s,
 | string scalar mode0, string scalar del)
{
	string scalar line, mode, m
	real scalar i, j, fh

	if (args()<4) del = char(9)
	mode = mm_strexpand(mode0, ("append", "replace"))
	m = "w"
	if (mode=="replace") unlink(fn)
	else if (mode=="append") m = "a"
	fh = fopen(fn, m)
	for (i=1; i<=rows(s); i++) {
		line = J(1,1,"")
		for (j=1; j<=cols(s); j++) {
			line = line + s[i,j]
			if (j<cols(s)) line = line + del
		}
		fput(fh, line)
	}
	fclose(fh)
}
end
