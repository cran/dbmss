\name{mhat}
\alias{mhat}
\title{
  Estimation of the m function
}
\description{
  Estimates the \emph{m} function
}
\usage{
mhat(X, r = NULL, ReferenceType, NeighborType = ReferenceType,
    CaseControl = FALSE, adjust = 1, CheckArguments = TRUE)
}
\arguments{
  \item{X}{
    A weighted, marked planar point pattern (\code{\link{wmppp.object}}).
  }
  \item{r}{
    A vector of distances. If \code{NULL}, a default value is set: 512 equally spaced values are used up to the median distance between points (following Duranton and Overman, 2005).
  }
  \item{ReferenceType}{
    One of the point types.
  }
  \item{NeighborType}{
    One of the point types. By default, the same as reference type.
  }
  \item{CaseControl}{
    Logical; if \code{TRUE}, the case-control version of \emph{M} is computed. \emph{ReferenceType} points are cases, \emph{NeighborType} points are controls.
  }
  \item{adjust}{
    Force the automatically selected bandwidth to be multiplied by \code{adjust}. Setting it to values lower than one (1/2 for example) will sharpen the estimation.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time in simulations for example, when the arguments have been checked elsewhere.
  }
}
\details{
  \emph{m} is a weighted, density, relative measure of a point pattern structure. Its value at any distance is the ratio of neighbors of the \emph{NeighborType} to all points around \emph{ReferenceType} points, normalized by its value over the windows.
  
  The number of neighbors at each distance is estimated by a Gaussian kernel whose bandwith is chosen optimally according to Sheather and Jones (1991). It can be sharpened or smoothed by multiplying it by \code{adjust}.
}
\value{
  An object of class \code{fv}, see \code{\link{fv.object}}, which can be plotted directly using \code{\link{plot.fv}}.
}
\references{
  Duranton, G. and Overman, H. G. (2005). Testing for Localisation Using Micro-Geographic Data. \emph{Review of Economic Studies} 72(4): 1077-1106.

  Sheather, S. J. and Jones M. C. (1991) A reliable data-based bandwidth selection method for kernel density estimation. Journal of the Royal Statistical Society B, 53(3): 683-690.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{mEnvelope}}, \code{\link{Kdhat}}
}
\examples{

data(paracou16)
plot(paracou16)

# Calculate M
plot(mhat(paracou16, , "V. Americana", "Q. Rosea"))
}
