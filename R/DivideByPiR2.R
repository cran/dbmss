DivideByPiR2 <-
function (Value, r) {
  NormalizedValue <- Value/pi/r^2
  NormalizedValue[r==0] <- NA              
  return(NormalizedValue)
}
