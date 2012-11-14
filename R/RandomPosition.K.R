RandomPosition.K <-
function(X) {
  RandomizedX <- runifpoint(X$n, win=X$window)
  marks(RandomizedX) <- data.frame(PointWeight = X$marks$PointWeight, PointType = X$marks$PointType)
  return(RandomizedX)
}
