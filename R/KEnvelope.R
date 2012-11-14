KEnvelope <-
function(NumberOfSimulations, Alpha, X, r, ReferenceType="", NeighborType="", SimulationType="RandomPosition") {
  # Compute simulations
  KSims <- t(replicate(NumberOfSimulations, SimulateK(X, r, ReferenceType, NeighborType, SimulationType)))
  # Compute the confidence envelope
  Envelope <- sapply(1:length(r), CriticalValues, KSims, Alpha)
  # Return the simulations and the envelope
  return(list(Simulations=KSims, Min=Envelope[1, ], Max=Envelope[2, ]))
}
