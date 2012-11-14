Lmm.r <-
function(X, r, ReferenceType="") {
  return(KtoL(Kmm.r(X, r, ReferenceType), r))
}
