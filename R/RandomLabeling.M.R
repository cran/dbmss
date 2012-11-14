RandomLabeling.M <-
function(X) {
  RandomizedX <- rlabel(X)
  RandomizedX$marks$PointWeight <- X$marks$PointWeight
  return(RandomizedX)
}
