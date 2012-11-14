PopulationIndependence.M <-
function(X, ReferenceType) {
  ReferencePP <- X[X$marks$PointType==ReferenceType]
  OtherPointsPP <- X[X$marks$PointType!=ReferenceType]
  RandomizedX <- superimpose(ReferencePP, rlabel(OtherPointsPP))
  return(RandomizedX)
}
