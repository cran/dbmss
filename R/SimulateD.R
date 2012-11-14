SimulateD <-
function(X, r, Cases, Controls, Intertype=FALSE) {
  # The only null hypothesis is random labelling (equivalently, random location)
  SimulatedPP <- rlabel(X)
  return(D.r(SimulatedPP, r, Cases, Controls, Intertype))
}
