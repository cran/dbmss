\name{D.r}
\alias{D.r}
\title{
  Estimation of the \emph{D} function
}
\description{
	Estimates the \emph{D} function}
\usage{
D.r(X, r, Cases, Controls, Intertype = FALSE)
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
  \item{Cases}{
	One of the point types.
  }
  \item{Controls}{
	One of the point types.
  }
  \item{Intertype}{
	Logical; if \code{TRUE}, \emph{D} is computed as \emph{Di} in Marcon and Puech (2012). 
  }
}
\details{
	The \emph{Di} function allows comparing the structure of the cases to that of the controls around cases, that is to say the comparison is made around the same points. This has been advocated by Arbia et al. (2008) and formalized by Marcon and Puech (2012).
}
\value{
	A vector containing \emph{D(r)} values.
}
\references{
	Arbia, G., Espa, G. and Quah, D. (2008). A class of spatial econometric methods in the empirical analysis of clusters of firms in the space. \emph{Empirical Economics} 34(1): 81-103.\cr
	Diggle, P. J. and Chetwynd, A. G. (1991). Second-Order Analysis of Spatial Clustering for Inhomogeneous Populations. \emph{Biometrics} 47(3): 1155-1163.\cr
	Marcon, E. and F. Puech (2012). A typology of distance-based measures of spatial concentration. \emph{HAL SHS}. 00679993.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\note{
	The computation of \code{D.r} relies on spatstat functions \code{\link{Kest}} and \code{\link{Kcross}}.
}

\seealso{
	\code{\link{K.r}}, \code{\link{DEnvelope}}, \code{\link{Kest}}, \code{\link{Kcross}}
}
\examples{

data(paracou16)
plot(paracou16)

# Calculate D
r <- 0:30
ActualValues.X <- D.r(paracou16, r, "V. Americana", "Q. Rosea", Intertype = TRUE)

# Plot
PlotResults(r, DivideByPiR2(ActualValues.X, r), , ylab="D / (pi R^2)", ReferenceValue=0)

}
