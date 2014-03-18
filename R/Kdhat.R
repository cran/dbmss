Kdhat <-
function(X, r = NULL, ReferenceType, NeighborType = ReferenceType, Weighted = FALSE, Original = TRUE, CheckArguments = TRUE) {
  
  if (CheckArguments) {
    CheckdbmssArguments()
  }
  
  # Select the bandwith: original choice by Duranton and Overman or optimized one.
  if (Original) {
    bw <- "nrd0"
  } else {
    bw <- "sj"
  }
  
  # Vectors to recognize point types
  IsReferenceType <- X$marks$PointType==ReferenceType
  IsNeighborType <- X$marks$PointType==NeighborType
  # Eliminate useless points
  Y <- X[IsReferenceType | IsNeighborType]
  # Update for Y
  IsReferenceType <- Y$marks$PointType==ReferenceType
  IsNeighborType <- Y$marks$PointType==NeighborType
  
  # Prepare a vector for distances between all point pairs.
  if (ReferenceType == NeighborType) {
    # Univariate Kd: n(n-1)/2 pairs
    NbDist <- sum(IsReferenceType)*(sum(IsReferenceType)-1)/2
  } else {
    # Bivariate Kd: n1*n2/2 pairs
    NbDist <- sum(IsReferenceType)*sum(IsNeighborType)/2
  }
  Dist <- vector(mode="double", length=NbDist)
  
  # Prepare a vector for weights if Weighted. Else, set a single value.
  if (Weighted) {
    Weights <- vector(mode="double", length=NbDist)
  } else {
    Weights <- 1
  }
  
  # C++ routine to fill distances and weights
  DistKd(Y$x, Y$y, Y$marks$PointWeight, Weights, Dist, IsReferenceType, IsNeighborType)
  
  # Estimate probability density.
  if (Weighted) {
    Density <- density(Dist, weights=Weights/sum(Weights), cut=0, bw=bw)
  } else {
    Density <- density(Dist, cut=0, bw=bw)  
  }

  if(is.null(r)) {
    # Return estimated values up to half the max distance (following D&O, 2005)
    r <- Density$x[1:256]
    Kd <- Density$y[1:256]
  } else {
    # Interpolate results at the chosen R
    Kd <- approx(Density$x, Density$y, xout=r)$y    
  }
  KdEstimate <- data.frame(r, Kd)
  colnames(KdEstimate) <- c("r", "Kd")
  
  # Return the values of g(r)
  return (fv(KdEstimate, argu="r", ylab=quote(Kd(r)), valu="Kd", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "Estimated Kd(r)"), unitname=X$window$unit, fname="Kd"))
}
