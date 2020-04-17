*! version 1.0.2  Ben Jann  05jan2008
version 9.0
mata:

string scalar mm_invtokens(string vector In, | real scalar noclean)
{
    string scalar Out
    string scalar tok
    real scalar i

    if (args()<2) noclean = 0
    Out = ""
    for (i=1; i<=length(In); i++) {
        tok = In[i]
        if (noclean)                 tok = "`" + `"""' + tok + `"""' + "'"
        else if (strpos(tok, `"""')) tok = "`" + `"""' + tok + `"""' + "'"
        else if (strpos(tok, " "))   tok = `"""' + tok + `"""'
        else if (tok=="")            tok = `"""' + `"""'
        if (i>1) tok = " " + tok
        Out = Out + tok
    }
    return(Out)
}
end
