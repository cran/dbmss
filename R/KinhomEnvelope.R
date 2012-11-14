KinhomEnvelope <-
function(NumberOfSimulations, Alpha, X, r, ReferenceType="", SimulationType="RandomPosition", lambda=NULL) {
  # Estimate intensity if it has not been provided.
  if (is.null(lambda)) {
    if (ReferenceType == "") {
      X.reduced <- X
    } else {
      X.reduced <- X[X$marks$PointType==ReferenceType]
    }
    lambda <- density.ppp(X.reduced, sigma=bw.diggle(X.reduced))
  }
  # Compute simulations 
  KinhomSims <- t(replicate(NumberOfSimulations, SimulateKinhom(X, r, ReferenceType, SimulationType, lambda)))
  # Compute the confidence envelope
  Envelope <- sapply(1:length(r), CriticalValues, KinhomSims, Alpha)
  # Return the simulations and the envelope
  return(list(Simulations=KinhomSims, Min=Envelope[1, ], Max=Envelope[2, ]))
}
