\name{Simulateg}
\alias{Simulateg}
\title{
	Simulations of point patterns to obtain values of g under the null hypothesis
}
\description{
	This function is used internally by \code{\link{gEnvelope}} to fill the simulation matrix.
}
\usage{
	Simulateg(X, r, ReferenceType = "", NeighborType = "", SimulationType = "RandomPosition")
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
	One of the point types.
  }
  \item{NeighborType}{
	One of the point types.
  }
  \item{SimulationType}{
	A string describing the null hypothesis to simulate. The null hypothesis may be "RandomPosition": points are drawn in a Poisson process; "RandomLabeling": randomizes point types, keeping locations unchanged; "PopulationIndependence": keeps reference points unchanged, randomizes other point locations.
  }
}

\value{
	A vector containing simulated g(r) values.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
	\code{\link{gEnvelope}}
}
\keyword{ internal }
