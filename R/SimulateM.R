SimulateM <-
function(X, r, ReferenceType, NeighborType, SimulationType="RandomLocation", CaseControl=FALSE) {
  SimulatedPP <- switch (SimulationType,
    RandomLocation = rlabel(X),
    RandomLabeling = RandomLabeling.M(X),
    PopulationIndependence = PopulationIndependence.M(X, ReferenceType)
    )
  return(M.r(SimulatedPP, r, ReferenceType, NeighborType, CaseControl))
}
