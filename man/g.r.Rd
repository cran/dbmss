\name{g.r}
\alias{g.r}
\title{
  Estimation of the g function
}
\description{
	Estimates the g function
}
\usage{
g.r(X, r, ReferenceType = "", NeighborType = "")
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
	one of the point types. Default is all point types.
  }
  \item{NeighborType}{
	one of the point types. Default is all point types.
  }
}
\details{
	The computation of g.r() relies on spatstat function markcorrint.
}
\value{
	A vector containing K(r) values.
}
\references{
	Stoyan, D. and Stoyan, H. (1994) Fractals, random shapes and point fields: methods of geometrical statistics. John Wiley and Sons. 
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}

\seealso{
	\code{\link{gEnvelope}}, \code{\link{pcf}}
}
\examples{

data(paracou16)
plot(paracou16)

# Calculate g
r <- 0:40
ActualValues.X <- g.r(paracou16, r, "Q. Rosea", "V. Americana")

# Plot
PlotResults(r, ActualValues.X, ylab="g", ReferenceValue=1)

}
