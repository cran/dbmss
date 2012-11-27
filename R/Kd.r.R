Kd.r <-
function(X, r, ReferenceType, NeighborType, Weighted=FALSE, Original=TRUE) {
  # Select the bandwith: original choice by Duranton and Overman or optimized one.
  if (Original) {
    bw = "nrd0"
  } else {
    bw = "sj"
  }
  
  # Compute the matrix of distances
  Dist <- pairdist.ppp(X)
  # Eliminate self point pair
  diag(Dist) <- NA
  
  # Vectors to recognize point types
  IsReferenceType <- X$marks$PointType==ReferenceType
  IsNeighborType <- X$marks$PointType==NeighborType

  # Reduce the matrix to pairs of interest
  Dist <- Dist[IsReferenceType, IsNeighborType]
  
  if (Weighted) {
    Weights <- matrix(rep(X$marks$PointWeight, each=X$n), nrow=X$n)
    Weights <- Weights[IsReferenceType, IsNeighborType]
    Density <- density(Dist, weights=Weights/sum(Weights), from=0, na.rm=TRUE, bw = bw)
  } else {
    Density <- density(Dist, from=0, na.rm=TRUE, bw = bw)
  }  
  # Interpolate results at the chosen R
  return(approx(Density$x, Density$y, xout=r)$y)
}
