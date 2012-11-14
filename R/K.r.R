K.r <-
function(X, r, ReferenceType="", NeighborType="") {
  # Eliminate erroneous configurations
  if ((ReferenceType=="" | NeighborType=="") & (ReferenceType!=NeighborType)) {
    stop("Either two or no point type must be specified.")
  }
  
  # K intratype calls Kest with the best edge-effect correction and returns the values
  Kintra <- function (X, r) {
    Kest.X <- Kest(X, r=r, correction="best")
    if (is.null(Kest.X$iso)) {
      return(Kest.X$trans)
    } else {
       return(Kest.X$iso)
    }
  } 
  
  # K intra
  if (ReferenceType=="" & NeighborType=="") {
    return(Kintra(X, r))
  }
  # K intra for a single point type
  if (ReferenceType==NeighborType) {
    X.reduced <- X[X$marks$PointType==ReferenceType]
    return(Kintra(X.reduced, r))
  }  
  # K inter calls Kcross. The marks must contain the type, with no weight.
  if (ReferenceType!=NeighborType) {
    X.cross <- X
    X.cross$marks <- X$marks$PointType
    Kcross.X <- Kcross(X.cross, i=ReferenceType , j=NeighborType, r=r, correction="best")
    if (is.null(Kcross.X$iso)) {
      return(Kcross.X$trans)
    } else {
      return(Kcross.X$iso)
    }

  }  
}
