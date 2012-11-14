D.r <-
function(X, r, Cases, Controls, Intertype=FALSE) {
  KCases <- K.r(X, r, Cases, Cases)            
  if (Intertype) {
    KControls <- K.r(X, r, Cases, Controls)   
  } else {
    KControls <- K.r(X, r, Controls, Controls)     
  }
  return(KCases-KControls)
}
