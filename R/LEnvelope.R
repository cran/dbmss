LEnvelope <-
function(NumberOfSimulations, Alpha, X, r, ReferenceType="", NeighborType="", SimulationType="RandomPosition") {
  Envelope <- KEnvelope(NumberOfSimulations, Alpha, X, r, ReferenceType, NeighborType, SimulationType)
  return(lapply(Envelope, KtoL, r))
}
