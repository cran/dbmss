GoFtest <- function(Envelope) {

  # Verify Envelope
  if (!inherits(Envelope, "envelope")) {
    stop("Envelope is not of class envelope")
  }
  # Verify simulations
  if (is.null(attr(Envelope, "simfuns"))) {
    stop("Envelope does not contain simulations in its attribute simfuns")
  } else {
    r <- as.data.frame(attr(Envelope, "simfuns"))[, 1]
    ActualValues <- Envelope$obs
    SimulatedValues <- as.data.frame(attr(Envelope, "simfuns"))[, -1]
  }

  NumberOfSimulations <- dim(SimulatedValues)[2]
  AverageSimulatedValues <- apply(SimulatedValues, 1, sum) / (NumberOfSimulations - 1)
  rIncrements <- (r - c(0,r)[seq_along(r)])[-1]

  # Ui calculate the statistic for a simulation
  Ui <- function(SimulationNumber) {
    Departure <- (SimulatedValues[, SimulationNumber] -
      AverageSimulatedValues)[seq_along(r) - 1]
    WeightedDeparture <- (Departure[!is.nan(Departure)])^2 *
      rIncrements[!is.nan(Departure)]
    return(sum(WeightedDeparture))
  }

  # Calculate the Ui statistic for all simulations
  SimulatedU <- vapply(
    seq_len(NumberOfSimulations),
    FUN = Ui,
    FUN.VALUE = 0
  )

  # Calculate the statistic for the actual value
  RecenteredValues <- (ActualValues - AverageSimulatedValues)[seq_along(r) - 1]
  WeightedRecenteredValues <- (RecenteredValues[!is.nan(RecenteredValues)])^2 *
    rIncrements[!is.nan(RecenteredValues)]
  ActualU <- sum(WeightedRecenteredValues)

  # Return the rank
  return(mean(ActualU < SimulatedU))
}
