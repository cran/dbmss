Kdhat <-
function(X, r = NULL, ReferenceType, NeighborType = ReferenceType, Weighted = FALSE,
         Original = TRUE, Approximate = ifelse(X$n < 10000, 0, 1), 
         CheckArguments = TRUE) {
  
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
  if (ReferenceType == "") {
    # All points (reference value as the center of the confidence interval)
    IsReferenceType <- IsNeighborType <- rep(TRUE, X$n)
    Y <- X
  } else {
    # Current use
    IsReferenceType <- X$marks$PointType==ReferenceType
    IsNeighborType <- X$marks$PointType==NeighborType    
    # Eliminate useless points
    Y <- X[IsReferenceType | IsNeighborType]
    # Update for Y
    IsReferenceType <- Y$marks$PointType==ReferenceType
    IsNeighborType <- Y$marks$PointType==NeighborType
  }

  if (Approximate) {
    # Round distances to save memory
    # Roughly estimate max distances
    if(is.null(r)) {
      # Default interval for R: between the min distance and the median one.
      rmin <- 0
      rmax <- sqrt((max(Y$x)-min(Y$x))^2 + (max(Y$y)-min(Y$y))^2)
    } else {
      rmin <- 0
      rmax <- max(r)
    }
    # Prepare 1024*Approximate steps
    rseq <- seq(from = rmin, to = rmax, length.out = 1024*Approximate)
    Nr <- length(rseq)
    # Prepare a matrix, single line, one value for each distance.
    NeighborWeights <- matrix(0.0, nrow=1, ncol=Nr)
    # Weights
    if (Weighted) {
      Weights <- Y$marks$PointWeight
    } else {
      Weights <- rep(1, Y$n)
    }
    
    # Call C routine to fill NeighborWeights
    CountNbdKd(rseq, Y$x, Y$y, Weights, NeighborWeights, IsReferenceType, IsNeighborType)
    
    # Adjust distances: values are the centers of intervals
    rseq <- c(0, (rseq[2:Nr]+rseq[1:Nr-1])/2)

    if(is.null(r)) {
      # Find the max distance, take the median
      rmax <- max(rseq[NeighborWeights>0])/2
    }
    
    # Estimate the density
    Density <- density(rseq, weights=NeighborWeights/sum(NeighborWeights), cut=0, to=rmax, bw=bw)  
    
  } else {
    # Classical estimation
    # Prepare a vector for distances between all point pairs.
    if (ReferenceType == NeighborType) {
      # Univariate Kd: n(n-1)/2 pairs
      NbDist <- sum(IsReferenceType)*(sum(IsReferenceType)-1)/2
    } else {
      # Bivariate Kd: n1*n2 pairs
      NbDist <- sum(IsReferenceType)*sum(IsNeighborType)
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
    
    if(is.null(r)) {
      # Default interval for R: between the min distance and the median one.
      rmin <- min(Dist)
      rmax <- median(Dist)
    } else {
      rmin <- min(r)
      rmax <- max(r)
    }
    # Estimate probability density.
    if (Weighted) {
      Density <- density(Dist, weights=Weights/sum(Weights), cut=0, from=rmin, to=rmax, bw=bw)
    } else {
      Density <- density(Dist, cut=0, from=rmin, to=rmax, bw=bw)  
    }
  }
  

  if(is.null(r)) {
    # Return estimated values
    r <- Density$x
    Kd <- Density$y
  } else {
    # Interpolate results at the chosen R
    Kd <- approx(Density$x, Density$y, xout=r)$y    
  }
  KdEstimate <- data.frame(r, Kd)
  colnames(KdEstimate) <- c("r", "Kd")
  
  # Return the values of Kd(r)
  return (fv(KdEstimate, argu="r", ylab=quote(Kd(r)), valu="Kd", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "Estimated Kd(r)"), unitname=X$window$unit, fname="Kd"))
}
