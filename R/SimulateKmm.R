SimulateKmm <-
function(X, r, ReferenceType="") {
  if (ReferenceType != "") {
    X.reduced <- X[X$marks$PointType==ReferenceType]
    SimulatedPP <- rlabel(X.reduced)
  } else {
    SimulatedPP <- rlabel(X)
  }
  return(Kmm.r(SimulatedPP, r))
}
