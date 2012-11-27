\name{PlotResults}
\alias{PlotResults}
\title{
	Plot of a distance-based measure against distances, with confidence envelopes
}
\description{
Provides a plot of a function against distance, with confidence interval envelopes.
}
\usage{
PlotResults(r, ActualValues, LocalCI, GlobalCI = NA, xlab = "r", ylab = "", 
    ReferenceValue = NA, Legend=FALSE, LegendItems=c("Value", "Local CI", "Global CI"), 
    LegendPosition="topright")
}
\arguments{
  \item{r}{
	A vector of distances.
  }
  \item{ActualValues}{
	A vector containing the values of a function (\emph{K}, \emph{M},...) at several distances.
  }
  \item{LocalCI}{
  A list (usually obtained by \code{\link{KEnvelope}} or similar) containing:\cr
	\code{Min}: A vector: the lower bound of the envelope\cr
	\code{Max}: A vector: the upper bound of the envelope
  }
  \item{GlobalCI}{
  A list (usually obtained by \code{\link{GlobalEnvelope}}) containing:\cr
	\code{Min}: A vector: the lower bound of the envelope\cr
	\code{Max}: A vector: the upper bound of the envelope
  }
  \item{xlab}{
  The label of the X axis of the plot
  }
  \item{ylab}{
  The label of the Y axis of the plot
  }
  \item{ReferenceValue}{
  The reference value of the function. A horizontal line representing the value is drawn.
  }
  \item{Legend}{
  Logical; if \code{TRUE}, a legend is added to the plot.
  }
  \item{LegendItems}{
  A vector of strings. Describes the function and its confidence envelopes.
  }
  \item{LegendPosition}{
  A single keyword from the list "\emph{bottomright}", "\emph{bottom}", "\emph{bottomleft}", "\emph{left}", "\emph{topleft}", "\emph{top}", "\emph{topright}", "\emph{right}" and "\emph{center}".
  }
}
\value{
  \code{NULL}
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\examples{

data(paracou16)
# Keep only 50% of points to run this example
X <- rthin(paracou16, 0.5)
plot(X)

# Calculate M
r <- seq(0, 30, 2)
ActualValues.X <- M.r(X, r, "V. Americana", "Q. Rosea")

# Calculate confidence envelope (should be 1000 simulations, reduced to 4 to save time)
NumberOfSimulations <- 4
Alpha <- .10
LocalEnvelope.X <- MEnvelope(NumberOfSimulations, Alpha, X, r, 
    "V. Americana", "Q. Rosea", "RandomLabeling")

# Plot
PlotResults(r, ActualValues.X, LocalEnvelope.X, ylab="M", ReferenceValue=1, 
    Legend=TRUE, LegendItems=c("M", "Local CI"))

}
