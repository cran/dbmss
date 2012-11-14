\name{Kd.r}
\alias{Kd.r}
\title{
  Estimation of the Kd function
}
\description{
	Estimates the Kd function
}
\usage{
Kd.r(X, r, ReferenceType, NeighborType, Weighted = FALSE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{Weighted}{
	if TRUE, estimates the Kemp function.
  }
}
\details{
	Kd is a density, relative measure of a point pattern structure.\cr
	Kd is computed efficiently by building a matrix of distances between point pairs and calculating the density of their distribution. As in Duranton and Overman (2005), the kernel estimator is Gaussian.
	The bandwdwith is calculated following Sheather and Jones (1991) rather than Silverman (1986: eq 3.31), see \code{\link{bw.SJ}} for more details.
}
\value{
	A vector containing M(r) values.
}
\references{
	Duranton, G. and Overman, H. G. (2005). Testing for Localisation Using Micro-Geographic Data. Review of Economic Studies 72(4): 1077-1106.\cr
	Sheather, S. J. and Jones, M. C. (1991) A reliable data-based bandwidth selection method for kernel density estimation. Journal of the Royal Statistical Society series B, 53, 683-690.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\note{
	Kd is computed using the density() function with a Gaussian kernel. This is more efficicent than Duranton and Overman's original algorithm.
}

\seealso{
	\code{\link{KdEnvelope}}, \code{\link{GlobalEnvelope}}, \code{\link{M.r}}
}
\examples{

data(paracou16)
plot(paracou16)

# Calculate Kd
r <- 0:100
ActualValues.X <- Kd.r(paracou16, r, "Q. Rosea", "Q. Rosea")

# Plot
PlotResults(r, ActualValues.X, ylab="Kd", Legend=TRUE, LegendItems="Kd", LegendPosition="bottomright")

}
