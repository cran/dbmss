mEnvelope <-
function(X, r = NULL, NumberOfSimulations = 100, Alpha = 0.05, ReferenceType, NeighborType = ReferenceType, 
         CaseControl = FALSE, adjust = 1, SimulationType = "RandomLocation", Global = FALSE) {
  
  CheckdbmssArguments()
  
  # Choose the null hypothesis
  SimulatedPP <- switch (SimulationType,
                         RandomLocation = expression(rRandomLocation(X, CheckArguments = FALSE)),
                         RandomLabeling = expression(rRandomLabelingM(X, CheckArguments = FALSE)),
                         PopulationIndependence = expression(rPopulationIndependenceM(X, ReferenceType, CheckArguments = FALSE))
  )
  if (is.null(SimulatedPP))
    stop(paste("The null hypothesis", sQuote(SimulationType), "has not been recognized."))
  # local envelope, keep extreme values for lo and hi (nrank=1)
  Envelope <- envelope(X, fun=mhat, nsim=NumberOfSimulations, nrank=1,
                       r=r, ReferenceType=ReferenceType, NeighborType=NeighborType, CaseControl=CaseControl, adjust=adjust,
                       CheckArguments = FALSE,
                       simulate=SimulatedPP, savefuns=TRUE
  )
  attr(Envelope, "einfo")$H0 <- switch (SimulationType,
                                        RandomLocation = "Random Location",
                                        RandomLabeling = "Random Labeling",
                                        PopulationIndependence = "Population Independence"
  )
  # Calculate confidence intervals
  Envelope <- FillEnveloppe(Envelope, Alpha, Global)
  # No edge effect correction
  attr(Envelope, "einfo")$valname <- NULL
  # Return the envelope
  return (Envelope)
}
