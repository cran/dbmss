\name{Mhat}
\alias{Mhat}
\title{
  Estimation of the M function
}
\description{
  Estimates the \emph{M} function
}
\usage{
Mhat(X, r, ReferenceType, NeighborType, CaseControl = FALSE, CheckArguments = TRUE)
}
\arguments{
  \item{X}{
    A weighted, marked planar point pattern (\code{\link{wmppp.object}}).
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
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time in simulations for example, when the arguments have been checked elsewhere.
  }
}
\details{
  \emph{M} is a weighted, cumulative, relative measure of a point pattern structure. Its value at any distance is the ratio of neighbors of the \emph{NeighborType} to all points around \emph{ReferenceType} points, normalized by its value over the windows.
}
\value{
  An object of class \code{fv}, see \code{\link{fv.object}}, which can be plotted directly using \code{\link{plot.fv}}.\cr
}
\references{
  Marcon, E. and Puech, F. (2010). Measures of the Geographic Concentration of Industries: Improving Distance-Based Methods. \emph{Journal of Economic Geography} 10(5): 745-762.\cr
  Marcon, E., F. Puech, et al. (2012). Characterizing the relative spatial structure of point patterns. \emph{International Journal of Ecology} 2012(Article ID 619281): 11.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{MEnvelope}}, \code{\link{Kdhat}}
}
\examples{

data(paracou16)
plot(paracou16)

# Calculate M
r <- seq(0, 30, 2)
plot(Mhat(paracou16, r, "V. Americana", "Q. Rosea"))
}