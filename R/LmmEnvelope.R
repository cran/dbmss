LmmEnvelope <-
function(NumberOfSimulations, Alpha, X, r, ReferenceType="") {
  Envelope <- KmmEnvelope(NumberOfSimulations, Alpha, X, r, ReferenceType)
  return(lapply(Envelope, KtoL, r))
}
