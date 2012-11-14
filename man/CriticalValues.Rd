\name{CriticalValues}
\alias{CriticalValues}
\title{
	Calcultate the critical values of a column of a matrix
}
\description{
	This function is used internally to apply quatiles() to each column of a matrix of simulations.
}
\usage{
CriticalValues(rNumber, Simulations, Alpha)
}
\arguments{
  \item{rNumber}{
	The number of the matrix column.
 }
  \item{Simulations}{
	The matrix of simulations.
 }
  \item{Alpha}{
	The risk level.
 }
}
\value{
	A vector of length 2 containg the Aplha/2 and 1-Alpha/2 quantiles.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\note{
	This function allows applying the quantile() function to a matrix.
}

\keyword{internal}
