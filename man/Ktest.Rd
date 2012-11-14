\name{Ktest}
\alias{Ktest}
\title{
  Test of a point pattern against Complete Spatial Randomness
}
\description{
  Tests the point pattern against CSR using values of the K function
}
\usage{
  Ktest(X, r)
}
\arguments{
  \item{X}{
  A point pattern (\code{\link{ppp.object}}). Marks are ignored. The window must be a rectangle sensu spatstat (tested by \code{\link{is.rectangle}}).
}
  \item{r}{
	A vector of distances.
}
}
\details{
	The test returns the risk to reject CSR erroneously, based on the distribution of the K function.
}
\value{
	A p-value.
}
\references{
	Lang, G. and Marcon, E. (2010). Testing randomness of spatial point patterns with the Ripley statistic. ArXiv e-prints 1006.1567.
}
\author{
	Gabriel Lang <Gabriel.Lang@agroparistech.fr>
}
\seealso{
	\code{\link{K.r}}, \code{\link{GoFtest}}
}
\examples{

# Simulate a Matern (Neyman Scott) point pattern
nclust <- function(x0, y0, radius, n) {
  return(runifdisc(n, radius, centre=c(x0, y0)))
}
X <- rNeymanScott(20, 0.1, nclust, radius=0.2, n=5)
plot(X)

# Test it
Ktest(X, r=seq(0.1, .5, .1))

}
