\name{Kd.r}
\alias{Kd.r}
\title{
  Estimation of the \emph{Kd} function
}
\description{
	Estimates the \emph{Kd} function
}
\usage{
Kd.r(X, r, ReferenceType, NeighborType, Weighted = FALSE, Original=TRUE)
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
	Logical; if \code{TRUE}, estimates the \emph{Kemp} function.
  }
  \item{Original}{
	Logical; if \code{TRUE} (by default), the original bandwidth selection by Duranton and Overman (2005) following Silverman (2006: eq 3.31) is used. If \code{FALSE}, it is calculated following Sheather and Jones (1991), ie. the state of the art. See \code{\link{bw.SJ}} for more details.
  }
}
\details{
	\emph{Kd} is a density, relative measure of a point pattern structure.\cr
	\emph{Kd} is computed efficiently by building a matrix of distances between point pairs and calculating the density of their distribution. The kernel estimator is Gaussian.
	The weighted \emph{Kd} function has been named \emph{Kemp} (\emph{emp} is for employees) by Duranton and Overman (2005).
}
\value{
	A vector containing \emph{Kd(r)} values.
}
\references{
	Duranton, G. and Overman, H. G. (2005). Testing for Localisation Using Micro-Geographic Data. \emph{Review of Economic Studies} 72(4): 1077-1106.\cr
	Sheather, S. J. and Jones, M. C. (1991) A reliable data-based bandwidth selection method for kernel density estimation. \emph{Journal of the Royal Statistical Society series B}, 53, 683-690.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
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
PlotResults(r, ActualValues.X, ylab="Kd", Legend=TRUE, LegendItems="Kd", 
    LegendPosition="bottomright")

}
