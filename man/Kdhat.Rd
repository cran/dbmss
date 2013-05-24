\name{Kdhat}
\alias{Kdhat}
\title{
  Estimation of the Kd function
}
\description{
  Estimates the \emph{Kd} function
}
\usage{
 Kdhat(X, r, ReferenceType, NeighborType = ReferenceType, Weighted = FALSE,
       Original = TRUE, CheckArguments = TRUE)
}
\arguments{
  \item{X}{
    A point pattern (\code{\link{wmppp.object}}).
  }
  \item{r}{
    A vector of distances. No default value is available.
  }
  \item{ReferenceType}{
    One of the point types.
  }
  \item{NeighborType}{
    One of the point types. By default, the same as reference type.
  }
  \item{Weighted}{
    Logical; if \code{TRUE}, estimates the \emph{Kemp} function.
  }
  \item{Original}{
    Logical; if \code{TRUE} (by default), the original bandwidth selection by Duranton and Overman (2005) following Silverman (2006: eq 3.31) is used. If \code{FALSE}, it is calculated following Sheather and Jones (1991), \emph{i.e.} the state of the art. See \code{\link{bw.SJ}} for more details.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time in simulations for example, when the arguments have been checked elsewhere.
  }
}
\details{
  \emph{Kd} is a density, absolute measure of a point pattern structure.\cr
  \emph{Kd} is computed efficiently by building a matrix of distances between point pairs and calculating the density of their distribution. The kernel estimator is Gaussian.
  The weighted \emph{Kd} function has been named \emph{Kemp} (\emph{emp} is for employees) by Duranton and Overman (2005).
}
\value{
  An object of class \code{fv}, see \code{\link{fv.object}}, which can be plotted directly using \code{\link{plot.fv}}.\cr
}
\references{
  Duranton, G. and Overman, H. G. (2005). Testing for Localisation Using Micro-Geographic Data. \emph{Review of Economic Studies} 72(4): 1077-1106.\cr
  Sheather, S. J. and Jones, M. C. (1991) A reliable data-based bandwidth selection method for kernel density estimation. \emph{Journal of the Royal Statistical Society series B}, 53, 683-690.\cr
  Silverman, B. W. (1986). \emph{Density estimation for statistics and data analysis}. Chapman and Hall, London.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{KdEnvelope}}, \code{\link{Mhat}}
}
\examples{

data(paracou16)
plot(paracou16)

# Calculate Kd
r <- 0:100
(Paracou <- Kdhat(paracou16, r, "Q. Rosea", "V. Americana"))
# Plot
plot(Paracou)
}
