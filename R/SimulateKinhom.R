SimulateKinhom <-
function(X, r, ReferenceType, SimulationType="RandomPosition", lambda) {
  if (SimulationType == "RandomLocation") {
    SimulatedPP <- rlabel(X)
    return(Kinhom.r(SimulatedPP, r, ReferenceType))  
  }
  if (SimulationType == "RandomPosition") {
    SimulatedPP <- rpoispp(lambda)
    # Points are simulated with no mark, Kinhom.r is called for all points.
    return(Kinhom.r(SimulatedPP, r, ReferenceType=""))    
  }
}
