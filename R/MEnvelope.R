MEnvelope <-
function(NumberOfSimulations, Alpha = .05, X, r, ReferenceType, NeighborType, SimulationType = "RandomLocation", CaseControl=FALSE) {
  # Warning for erroneous configurations
  if (CaseControl & (SimulationType != "RandomLocation")) {
    warning(paste("The null hypothesis for case-control M should be RandomLocation. The argument used is", SimulationType))
  }
  # Compute simulations
  Msims <- t(replicate(NumberOfSimulations, SimulateM(X, r, ReferenceType, NeighborType, SimulationType, CaseControl)))
  # Compute the local confidence envelope
  Envelope <- sapply(1:length(r), CriticalValues, Msims, Alpha)
  # Return the simulations and the envelope
  return(list(Simulations=Msims, Min=Envelope[1, ], Max=Envelope[2, ]))
}
