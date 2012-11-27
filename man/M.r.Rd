\name{M.r}
\alias{M.r}
\title{
  Estimation of the \emph{M} function
}
\description{
  Estimates the \emph{M} function
}
\usage{
  M.r(X, r, ReferenceType, NeighborType, CaseControl = FALSE)
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
  \item{CaseControl}{
	Logical; if \code{TRUE}, the case-control version of \emph{M} is computed. \emph{ReferenceType} points are cases, \emph{NeighborType} points are controls.
  }
}
\details{
	\emph{M} is a weighted, cumulative, relative measure of a point pattern structure. Its value at any distance is the ratio of neighbors of the \emph{NeighborType} to all points around \emph{ReferenceType} points, normalized by its value over the windows.
}
\value{
	A vector containing \emph{M(r)} values.
}
\references{
	Marcon, E. and Puech, F. (2010). Measures of the Geographic Concentration of Industries: Improving Distance-Based Methods. \emph{Journal of Economic Geography} 10(5): 745-762.\cr
	Marcon, E., F. Puech, et al. (2012). Characterizing the relative spatial structure of point patterns. \emph{International Journal of Ecology} 2012(Article ID 619281): 11.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
	\code{\link{MEnvelope}}, \code{\link{Kd.r}}
}
\examples{

data(paracou16)
plot(paracou16)

# Calculate M
r <- seq(0, 30, 2)
ActualValues.X <- M.r(paracou16, r, "V. Americana", "Q. Rosea")

# Plot
PlotResults(r, ActualValues.X, ylab="M", ReferenceValue=1)

}
