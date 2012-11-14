\name{SimulateKmm}
\alias{SimulateKmm}
\title{
	Simulations of point patterns to obtain values of Kmm under the null hypothesis
}
\description{
	This function is used internally by \code{\link{KmmEnvelope}} to fill the simulation matrix.
}
\usage{
SimulateKmm(X, r, ReferenceType = "")
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
  \item{ReferenceType}{
	One of the point types. Others are ignored. Default is all point types.
  }
}

\value{
	A vector containing simulated Kmm(r) values.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
	\code{\link{KmmEnvelope}}
}
\keyword{ internal }
