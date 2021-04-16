## ----global_options, include=FALSE--------------------------------------------
set.seed(2018)

## ----wmppp, warning=FALSE, message=FALSE--------------------------------------
library("dbmss")
# Draw the coordinates of 10 points
X <- runif(10)
Y <- runif(10)
# Draw the point types.
PointType   <- sample(c("A", "B"), 10, replace=TRUE)
# Plot the point pattern. Weights are set to 1 ant the window is adjusted
autoplot(wmppp(data.frame(X, Y, PointType)))

## ----paracou------------------------------------------------------------------
# Plot (second column of marks is Point Types) 
autoplot(paracou16, 
  labelSize = expression("Basal area (" ~cm^2~ ")"), 
  labelColor = "Species")

## ----m------------------------------------------------------------------------
autoplot(Mhat(paracou16, , "V. Americana", "Q. Rosea"), main="")

## -----------------------------------------------------------------------------
autoplot(KdEnvelope(paracou16, , ReferenceType="Q. Rosea", Global=TRUE), main="")

