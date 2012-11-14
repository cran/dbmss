DEnvelope <-
function(NumberOfSimulations, Alpha, X, r, Cases, Controls, Intertype=FALSE) {
  # Compute simulations
  DSims <- t(replicate(NumberOfSimulations, SimulateD(X, r, Cases, Controls, Intertype)))
  # Compute the confidence envelope
  Envelope <- sapply(1:length(r), CriticalValues, DSims, Alpha)
  # Return the simulations and the envelope
  return(list(Simulations=DSims, Min=Envelope[1, ], Max=Envelope[2, ]))  
}
