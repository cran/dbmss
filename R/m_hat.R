mhat <-
function(X, r = NULL, ReferenceType, NeighborType = ReferenceType, CaseControl = FALSE, adjust = 1, CheckArguments = TRUE) {
  # Eliminate erroneous configurations
  if (CheckArguments) {
    CheckdbmssArguments()
    if (CaseControl & (ReferenceType==NeighborType)) {
      warning("Cases and controls are identical.")
      return(rep(1,length(r)))
    }
  }
  
  # Vectors to recognize point types
  IsReferenceType <- X$marks$PointType==ReferenceType
  IsNeighborType <- X$marks$PointType==NeighborType
  
  # Global ratio
  if (ReferenceType==NeighborType | CaseControl) {
    WrMinusReferencePoint <- sum(X$marks$PointWeight[IsReferenceType])-X$marks$PointWeight
    Wn <- WrMinusReferencePoint[IsReferenceType]
  } else {
    Wn <- sum(X$marks$PointWeight[IsNeighborType])
  }
  if (CaseControl) {
    Wa <- sum(X$marks$PointWeight[IsNeighborType]) 
  } else {
    WaMinusReferencePoint <- sum(X$marks$PointWeight)-X$marks$PointWeight
    Wa <- WaMinusReferencePoint[IsReferenceType]
  }
  GlobalRatio <- Wn/Wa
  
  # Distance matrix
  pdAll <- pairdist(X)
  diag(pdAll) <- NA    
  # Choose the bandwith based on all distance pairs
  h <- bw.nrd0(pdAll[upper.tri(pdAll)])*adjust
  # Keep the lines of the matrix corresponding to reference points (cases).
  pdAll <- pdAll[IsReferenceType, ]
  pdNeighbors <- pdAll[, IsNeighborType]
  
  if(is.null(r)) {
    # Default interval for R: between the min distance and the median one.
    rmin <- min(pdAll, na.rm=TRUE)
    rmax <- median(pdAll, na.rm=TRUE)
  } else {
    rmin <- min(r)
    rmax <- max(r)
  }
  
  # Calculate densities of neighbors (with unnormalized weights so suppress warnings)
  Djc <- t(apply(pdNeighbors, 1, function(x) suppressWarnings(density(x, bw=h, weights=X$marks$PointWeight[IsNeighborType][!is.na(x)], cut=0, from=rmin, to=rmax, na.rm=TRUE))$y))
  Dj <- t(apply(pdAll, 1, function(x) suppressWarnings(density(x, bw=h, weights=X$marks$PointWeight[!is.na(x)], cut=0, from=rmin, to=rmax, na.rm=TRUE))$y))
  
  # Calculate the local ratio (at distance r)
  LocalRatio <- Djc/Dj
  # Divide it by the global ratio. Ignore points with no neighbor at all.
  Mvalues <- colSums(LocalRatio)/sum(GlobalRatio)
  
  # Interpolate if necessary
  if (is.null(r)) {
    r <- density(pdNeighbors[1,], bw=h, cut=0, from=rmin, to=rmax, na.rm=TRUE)$x
  } else {
    Mvalues <- approx(r, Mvalues, xout=r)$y 
  }
  # Put the results into an fv object
  MEstimate <- data.frame(r, rep(1, length(r)), Mvalues)
  colnames(MEstimate) <- c("r", "theo", "m")
  
  # Return the values of M(r)
  return (fv(MEstimate, argu="r", ylab=quote(m(r)), valu="m", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "%s[ind](r)", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "theoretical independent m(r)", "Estimated m(r)"), unitname=X$window$unit, fname="m")) 
}
