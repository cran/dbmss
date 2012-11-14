Simulateg <-
function(X, r, ReferenceType="", NeighborType="", SimulationType="RandomPosition") {
  SimulatedPP <- switch (SimulationType,
    RandomPosition = RandomPosition.K(X),
    RandomLabeling = rlabel(X),
    PopulationIndependence = PopulationIndependence.K(X, ReferenceType, NeighborType)
    )
  return(g.r(SimulatedPP, r, ReferenceType, NeighborType))
}
