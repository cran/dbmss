Dhat <-
function(X, r = NULL, Cases, Controls, Intertype = FALSE, CheckArguments = TRUE) {
  if (CheckArguments) {
    CheckdbmssArguments()
  }
  # K of cases.
  KCases <- Khat(X, r, Cases, Cases, CheckArguments)   
  # K of controls. r must be those of cases.
  if (Intertype) {
    KControls <- Khat(X, KCases$r, Cases, Controls, CheckArguments)   
  } else {
    KControls <- Khat(X, KCases$r, Controls, Controls, CheckArguments)     
  }
  # Calculate the difference (a difference between fv's yields a dataframe)
  Dvalues <- KCases-KControls
  
  return (fv(cbind(as.data.frame(KCases)[1], Dvalues[2:3]), argu="r", ylab=quote(D(r)), valu=attr(KCases, "valu"), fmla=attr(KCases, "fmla"), alim=attr(KCases, "alim"), labl=c("r", "%s[ind](r)", "hat(%s)[iso](r)"), desc=attr(KCases, "desc"), unitname=attr(KCases, "unitname"), fname="D"))
}
