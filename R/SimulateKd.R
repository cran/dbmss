SimulateKd <-
function(X, r, ReferenceType, NeighborType, Weighted=FALSE, SimulationType="RandomLocation") {
  SimulatedPP <- switch (SimulationType,
    RandomLocation = rlabel(X),
    RandomLabeling = RandomLabeling.M(X),
    PopulationIndependence = PopulationIndependence.M(X, ReferenceType)
    )
  return(Kd.r(SimulatedPP, r, ReferenceType, NeighborType, Weighted))
}
