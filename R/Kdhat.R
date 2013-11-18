Kdhat <-
function(X, r, ReferenceType, NeighborType = ReferenceType, Weighted = FALSE, Original = TRUE, CheckArguments = TRUE) {
  
  if (CheckArguments) {
    CheckdbmssArguments()
  }
  
  # Select the bandwith: original choice by Duranton and Overman or optimized one.
  if (Original) {
    bw = "nrd0"
  } else {
    bw = "sj"
  }
  
  # Vectors to recognize point types
  IsReferenceType <- X$marks$PointType==ReferenceType
  IsNeighborType <- X$marks$PointType==NeighborType
  # Eliminate useless points
  X <- X[IsReferenceType | IsNeighborType]
  # Compute the matrix of distances
  Dist <- pairdist.ppp(X)
  # Eliminate self point pair
  diag(Dist) <- NA
    
  # Reduce the matrix to pairs of interest
  IsReferenceType <- X$marks$PointType==ReferenceType
  IsNeighborType <- X$marks$PointType==NeighborType
  Dist <- Dist[IsReferenceType, IsNeighborType]
  
  if (Weighted) {
    Weights <- matrix(rep(X$marks$PointWeight, each=X$n), nrow=X$n)
    Weights <- Weights[IsReferenceType, IsNeighborType]
    # Eliminate self point pair so that density works later
    if (NeighborType == ReferenceType) {
      diag(Weights) <- NA
      Dist <- Dist[!is.na(Dist)]
      Weights <- Weights[!is.na(Weights)]
    }
    Density <- density(Dist, weights=Weights/sum(Weights), from=0, bw=bw)
  } else {
    Density <- density(Dist, from=0, na.rm=TRUE, bw=bw)
  }  
  # Interpolate results at the chosen R
  Kd <- approx(Density$x, Density$y, xout=r)$y
  KdEstimate <- data.frame(r, Kd)
  colnames(KdEstimate) <- c("r", "Kd")
  
  # Return the values of g(r)
  return (fv(KdEstimate, argu="r", ylab=quote(Kd(r)), valu="Kd", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "Estimated Kd(r)"), unitname=X$window$unit, fname="Kd"))
}
