Kinhom.r <-
function(X, r, ReferenceType="", lambda=NULL) {
  
  # K intratype calls Kest with the best edge-effect correction and returns the values
  Kiintra <- function (X, r, lambda) {
    # Estimate intensity if it has not been provided
    if (is.null(lambda)) {
      lambda <- as.vector(density.ppp(X, sigma=bw.diggle(X), at="points"))
    }
    # Calculate Kinhom according to lambda
    Kinhom.X <- Kinhom(X, r=r, correction=c("isotropic", "translate", "border"), normpower=2, lambda=lambda)
    # Return the best edge-effect corrected value : iso if possible, border is the worst
    if (is.null(Kinhom.X$iso)) {
      if (is.null(Kinhom.X$trans)) {
        return(Kinhom.X$border)
      } else {
        return(Kinhom.X$trans)
      }
    } else {
       return(Kinhom.X$iso)
    }
  } 
  
  # K intra
  if (ReferenceType=="") {
    return(Kiintra(X, r, lambda))
  }
  # K intra for a single point type
  if (ReferenceType!="") {
    X.reduced <- X[X$marks$PointType==ReferenceType]
    return(Kiintra(X.reduced, r, lambda))
  }  

}
