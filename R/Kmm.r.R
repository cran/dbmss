Kmm.r <-
function(X, r, ReferenceType="") {

  # KmmBymarkcorrint calls markcorrint with the best edge-effect correction and returns the values
  KmmBymarkcorrint <- function (X, r) {
    X.marked <- X
    # Weights are normalized so that their mean is 1 because markcorrint returns Kmm * mean weight instead of Kmm (as of v. 1.27-0 od spatstat).
    X.marked$marks <- X$marks$PointWeight/mean(X$marks$PointWeight)
    Kest.X <- markcorrint(X.marked, correction="best")
    if (is.null(Kest.X$iso)) {
      return(approx(Kest.X$r, Kest.X$trans, xout=r)$y)
    } else {
      return(approx(Kest.X$r, Kest.X$iso, xout=r)$y)
    }
  } 
  
  # Kmm all points or specified point type
  if (ReferenceType=="") {
    return(KmmBymarkcorrint(X, r))
    } else {
    X.reduced <- X[X$marks$PointType==ReferenceType]
    return(KmmBymarkcorrint(X.reduced, r))
  }   
}
