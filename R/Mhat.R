Mhat <-
function(X, r, ReferenceType, NeighborType, CaseControl = FALSE, CheckArguments = TRUE) {
  # Eliminate erroneous configurations
  if (CheckArguments) {
    CheckdbmssArguments()
    if (CaseControl & (ReferenceType==NeighborType)) {
      warning("Cases and controls are identical.")
      return(rep(1,length(r)))
    }
  }
  
  # Compute the matrix of distances (squared to save time)
  Dist <- pairdist.ppp(X, squared=TRUE)
  
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
  
  # M.distance calculates M at a single distance
  M.distance <- function(distance) {
    # Select point pairs less than r apart
    IsCloseEnough <- (Dist <= distance)
    # Eliminate point pairs made of a single point
    diag(IsCloseEnough) <- FALSE
    # Calculate the weight of neighbors
    IsCloseEnough <- IsCloseEnough[IsReferenceType, ]
    if (CaseControl) {
      IsCloseEnoughAndCase <- t(t(IsCloseEnough) & IsReferenceType)
      NeighborTypeWeight <- IsCloseEnoughAndCase %*% X$marks$PointWeight
      IsCloseEnoughAndControl <- t(t(IsCloseEnough) & IsNeighborType)
      AllNeighborWeight <- IsCloseEnoughAndControl %*% X$marks$PointWeight
    } else {                                                    
      IsCloseEnoughAndNeighborType <- t(t(IsCloseEnough) & IsNeighborType)
      NeighborTypeWeight <- IsCloseEnoughAndNeighborType %*% X$marks$PointWeight
      AllNeighborWeight <- IsCloseEnough %*% X$marks$PointWeight
    }
    # Calculate the local ratio
    LocalRatio <- NeighborTypeWeight/AllNeighborWeight
    # Calculate M, eliminate undefined values (no neighbor of any type)
    return(sum(LocalRatio[is.finite(LocalRatio)])/sum(GlobalRatio[is.finite(LocalRatio)]))
  }
    
  # Build a dataframe with r, theoretical value = 1, and M(r)
  MEstimate <- data.frame(r, rep(1, length(r)), sapply(r*r, M.distance))
  colnames(MEstimate) <- c("r", "theo", "M")
  
  # Return the values of M(r)
  return (fv(MEstimate, argu="r", ylab=quote(M(r)), valu="M", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "%s[ind](r)", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "theoretical independent M(r)", "Estimated M(r)"), unitname=X$window$unit, fname="M"))
  
}
