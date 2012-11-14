KmmEnvelope <-
function(NumberOfSimulations, Alpha, X, r, ReferenceType="") {
  # Compute simulations
  KmmSims <- t(replicate(NumberOfSimulations, SimulateKmm(X, r, ReferenceType)))
  # Compute the confidence envelope
  Envelope <- sapply(1:length(r), CriticalValues, KmmSims, Alpha)
  # Return the simulations and the envelope
  return(list(Simulations=KmmSims, Min=Envelope[1, ], Max=Envelope[2, ]))
}
