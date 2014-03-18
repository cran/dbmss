\name{KdEnvelope}
\alias{KdEnvelope}
\title{
  Estimation of the confidence envelope of the Kd function under its null hypothesis
}
\description{
    Simulates point patterns according to the null hypothesis and returns the envelope of \emph{Kd} according to the confidence level.
}
\usage{
KdEnvelope(X, r = NULL, NumberOfSimulations = 100, Alpha = 0.05, ReferenceType, 
           NeighborType = ReferenceType, Weighted = FALSE, Original = TRUE, 
           SimulationType = "RandomLocation", Global = FALSE)
}
\arguments{
  \item{X}{
    A point pattern (\code{\link{wmppp.object}}).
  }
  \item{r}{
    A vector of distances. If \code{NULL}, a default value is set: 512 equally spaced values are used, and the first 256 are returned, corresponding to half the maximum distance between points (following Duranton and Overman, 2005).
  }
  \item{NumberOfSimulations}{
    The number of simulations to run.
  }
  \item{Alpha}{
    The risk level.
  }
  \item{ReferenceType}{
    One of the point types. By default, the same as reference type.
  }
  \item{NeighborType}{
    One of the point types.
  }
  \item{Weighted}{
    Logical; if \code{TRUE}, estimates the \emph{Kemp} function.
  }
  \item{Original}{
    Logical; if \code{TRUE} (by default), the original bandwidth selection by Duranton and Overman (2005) following Silverman (2006: eq 3.31) is used. If \code{FALSE}, it is calculated following Sheather and Jones (1991), \emph{i.e.} the state of the art. See \code{\link{bw.SJ}} for more details.
  }
  \item{SimulationType}{
    A string describing the null hypothesis to simulate. The null hypothesis may be
    "\emph{RandomLocation}": points are redistributed on the actual locations;
    "\emph{RandomLabeling}": randomizes point types, keeping locations and weights unchanged;
    "\emph{PopulationIndependence}": keeps reference points unchanged, randomizes other point locations.
  }
  \item{Global}{
    Logical; if \code{TRUE}, a global envelope sensu Duranton and Overman (2005) is calculated. 
  }
}
\details{
  This envelope is local by default, that is to say it is computed separately at each distance. See Loosmore and Ford (2006) for a discussion.\cr
  The global enveloppe is calculated by iteration: the simulations reaching one of the upper or lower values at any distance are eliminated at each step. The process is repeated until \emph{Alpha / Number of simulations} simulations are dropped. The remaining upper and lower bounds at all distances constitute the global envelope. Interpolation is used if the exact ratio cannot be reached.
}
\value{
  An envelope object (\code{\link{envelope}}). There are methods for print and plot for this class.\cr
  The \code{fv} contains the observed value of the function, its average simulated value and the confidence enveloppe.
}
\references{
  Duranton, G. and Overman, H. G. (2005). Testing for Localisation Using Micro-Geographic Data. \emph{Review of Economic Studies} 72(4): 1077-1106.\cr
  Kenkel, N. C. (1988). Pattern of Self-Thinning in Jack Pine: Testing the Random Mortality Hypothesis. \emph{Ecology} 69(4): 1017-1024.\cr
  Loosmore, N. B. and Ford, E. D. (2006). Statistical inference using the G or K point pattern spatial statistics. \emph{Ecology} 87(8): 1925-1931.\cr
  Marcon, E. and F. Puech (2012). A typology of distance-based measures of spatial concentration. \emph{HAL SHS}. 00679993.\cr
  Silverman, B. W. (1986). \emph{Density estimation for statistics and data analysis}. Chapman and Hall, London.

}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Kdhat}}
}

\examples{

data(paracou16)
plot(paracou16[paracou16$marks$PointType=="Q. Rosea"])

# Calculate confidence envelope
plot(KdEnvelope(paracou16, , ReferenceType="Q. Rosea", Global=TRUE))
}