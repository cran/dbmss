gEnvelope <-
function(NumberOfSimulations, Alpha, X, r, ReferenceType="", NeighborType="", SimulationType = "RandomPosition") {
  # Compute simulations
  gSims <- t(replicate(NumberOfSimulations, Simulateg(X, r, ReferenceType, NeighborType, SimulationType)))
  # Compute the confidence envelope
  Envelope <- sapply(1:length(r), CriticalValues, gSims, Alpha)
  # Return the simulations and the envelope
  return(list(Simulations=gSims, Min=Envelope[1, ], Max=Envelope[2, ]))
}
