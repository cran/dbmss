## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(
  fig.width = 5,       # Larger figures (default is 3, only legend is visible)
  out.width = "100%"
)
set.seed(2018)

## ----wmppp, warning=FALSE, message=FALSE--------------------------------------
library("dbmss")
# Draw the coordinates of 10 points
X <- runif(10)
Y <- runif(10)
# Draw the point types.
PointType <- sample(c("A", "B"), size = 10, replace = TRUE)
# Plot the point pattern. Weights are set to 1 ant the window is adjusted
autoplot(wmppp(data.frame(X, Y, PointType)))

## ----paracou------------------------------------------------------------------
# Plot (second column of marks is Point Types) 
autoplot(
  paracou16, 
  labelSize = expression("Basal area (" ~cm^2~ ")"), 
  labelColor = "Species"
)

## ----m------------------------------------------------------------------------
autoplot(
  Mhat(
    paracou16, 
    ReferenceType = "V. Americana", 
    NeighborType = "Q. Rosea"
  ), 
  main = ""
)

## -----------------------------------------------------------------------------
autoplot(
  KdEnvelope(paracou16, ReferenceType = "Q. Rosea", Global = TRUE), 
  main = ""
)

## -----------------------------------------------------------------------------
# Calculate individual intertype M(distance) value
ReferenceType <- "V. Americana"
NeighborType <- "Q. Rosea"
fvind <- Mhat(
  paracou16, 
  r = c(0, 30), 
  ReferenceType = ReferenceType, 
  NeighborType = NeighborType, 
  Individual = TRUE
)
# Plot the point pattern with values of M(30 meters)
p16_map <- Smooth(
  paracou16, 
  fvind = fvind, 
  distance = 30,
  # Resolution
  Nbx = 512, 
  Nby = 512
)
par(mar = rep(0, 4))
plot(p16_map, main = "")
# Add the reference points to the plot
is.ReferenceType <- marks(paracou16)$PointType == ReferenceType
points(
  x = paracou16$x[is.ReferenceType], 
  y = paracou16$y[is.ReferenceType], 
  pch = 20
)
# Add contour lines
contour(p16_map, nlevels = 5, add = TRUE)

