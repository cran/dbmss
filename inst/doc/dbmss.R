## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(
  fig.width = 5       # Larger figures (default is 3, only legend is visible)
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

