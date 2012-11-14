KdEnvelope <-
function(NumberOfSimulations, Alpha, X, r, ReferenceType, NeighborType, Weighted=FALSE, SimulationType = "RandomLocation") {
  # Compute simulations
  KdSims <- t(replicate(NumberOfSimulations, SimulateKd(X, r, ReferenceType, NeighborType, Weighted, SimulationType)))
  # Compute the confidence envelope
  Envelope <- sapply(1:length(r), CriticalValues, KdSims, Alpha)
  # Return the simulations and the envelope
  return(list(Simulations=KdSims, Min=Envelope[1, ], Max=Envelope[2, ]))
}
