CriticalValues <-
function(rNumber, Simulations, Alpha) {
  return(quantile(Simulations[ , rNumber], probs=c(Alpha, 1-Alpha), na.rm=TRUE))
}
