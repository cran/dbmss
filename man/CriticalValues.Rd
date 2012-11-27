\name{CriticalValues}
\alias{CriticalValues}
\title{
	Calculate the critical values of a column of a matrix
}
\description{
	This function is used internally to apply \code{\link{quantile}} to each column of a matrix of simulations.
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
	A vector of length 2 containg the \eqn{Alpha/2} and \eqn{1-Alpha/2} quantiles.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\note{
	This function allows applying the \code{\link{quantile}} function to a matrix.
}

\keyword{internal}
