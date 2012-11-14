\name{SimulateD}
\alias{SimulateD}
\title{
	Simulations of point patterns to obtain values of D under the null hypothesis
}
\description{
	This function is used internally by \code{\link{DEnvelope}} to fill the simulation matrix.
}
\usage{
SimulateD(X, r, Cases, Controls, Intertype = FALSE)
}
\arguments{
  \item{X}{
  A point pattern (\code{\link{ppp.object}}), marks must be a dataframe with two columns:\cr
    \code{PointType}: labels, as factors.\cr
    \code{PointWeight}: weights.
  }
  \item{r}{
	A vector of distances.
  }
  \item{Cases}{
	One of the point types
  }
  \item{Controls}{
	One of the point types.
  }
  \item{Intertype}{
	if TRUE, D is computed as Di in Marcon and Puech (2012). 
  }
}
\value{
	A vector containing simulated D(r) values.
}
\references{
	Kenkel, N. C. (1988). Pattern of Self-Thinning in Jack Pine: Testing the Random Mortality Hypothesis. Ecology 69(4): 1017-1024.\cr
	Marcon, E. and F. Puech (2012). A typology of distance-based measures of spatial concentration. HAL SHS. 00679993.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
	\code{\link{DEnvelope}}
}
\keyword{ internal }
