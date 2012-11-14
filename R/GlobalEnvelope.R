GlobalEnvelope <-
function(Simulations, Alpha) {
  # Initialize
  NumberOfSimulations <- dim(Simulations)[1]
  TargetNumber <- (1-Alpha)*NumberOfSimulations
  KeptSimulations <- Simulations
  PresentMinValues <- apply(Simulations, 2, min)
  PresentMaxValues <- apply(Simulations, 2, max)
  # Loop until the target number of simulations is kept
  while(dim(KeptSimulations)[1] > TargetNumber) {
    # Remember previous min and max
    PreviousMinValues <- PresentMinValues
    PreviousMaxValues <- PresentMaxValues
    # Select the simulations that gave extreme values
    SimulationsToDrop <- c(unlist(apply(KeptSimulations, 2, which.min)), unlist(apply(KeptSimulations, 2, which.max)))
    # Drop them
    KeptSimulations <- KeptSimulations[-SimulationsToDrop, ]              
    # Calculate min and max
    PresentMinValues <- apply(KeptSimulations, 2, min)
    PresentMaxValues <- apply(KeptSimulations, 2, max)
  }  
  # Interpolate because the kept number of simulations is not always the target
  NumberOfKeptSimulations <- dim(KeptSimulations)[1]
  MinValues <- PresentMinValues + (PreviousMinValues-PresentMinValues)/NumberOfSimulations*(TargetNumber-NumberOfKeptSimulations)
  MaxValues <- PresentMaxValues + (PreviousMaxValues-PresentMaxValues)/NumberOfSimulations*(TargetNumber-NumberOfKeptSimulations)
  return(list(Min=MinValues, Max=MaxValues))
}
