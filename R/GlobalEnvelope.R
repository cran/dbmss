GlobalEnvelope <-
function(Simulations, Alpha) {
  
  verifyclass(Simulations, "fv")  
  
  # Initialize
  Simulations <- as.data.frame(Simulations)[, -1] # Eliminate r
  NumberOfSimulations <- dim(Simulations)[2]
  TargetNumber <- (1-Alpha)*NumberOfSimulations
  KeptSimulations <- Simulations
  PresentMinValues <- apply(Simulations, 1, min)
  PresentMaxValues <- apply(Simulations, 1, max)
  # Loop until the target number of simulations is kept
  while(dim(KeptSimulations)[2] > TargetNumber) {
    # Remember previous min and max
    PreviousMinValues <- PresentMinValues
    PreviousMaxValues <- PresentMaxValues
    # Select the simulations that gave extreme values
    SimulationsToDrop <- c(unlist(apply(KeptSimulations, 1, which.min)), unlist(apply(KeptSimulations, 1, which.max)))
    # Drop them
    KeptSimulations <- KeptSimulations[, -SimulationsToDrop]  
    # Fails if no simulations are left
    if (is.null(dim(KeptSimulations)))
      stop("Global envelope could not be calculated. More simulations are necessary.")
    # Calculate min and max
    PresentMinValues <- apply(KeptSimulations, 1, min)
    PresentMaxValues <- apply(KeptSimulations, 1, max)
  }  
  # Interpolate because the kept number of simulations is not always the target
  NumberOfKeptSimulations <- dim(KeptSimulations)[2]
  Glo <- PresentMinValues + (PreviousMinValues-PresentMinValues)/NumberOfSimulations*(TargetNumber-NumberOfKeptSimulations)
  Ghi <- PresentMaxValues + (PreviousMaxValues-PresentMaxValues)/NumberOfSimulations*(TargetNumber-NumberOfKeptSimulations)
  return(rbind(Glo, Ghi))
}
