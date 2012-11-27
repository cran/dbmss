\name{Kmm.r}
\alias{Kmm.r}
\title{
  Estimation of the \emph{Kmm} function
}
\description{
  Estimates of the \emph{Kmm} function
}
\usage{
Kmm.r(X, r, ReferenceType = "")
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
	One of the point types. Others are ignored. Default is all point types.
  }
}
\details{
	The \emph{Kmm} function is used to test the independence of marks.
}
\value{
	A vector containing \emph{Kmm(r)} values.
}
\references{
	Penttinen, A., Stoyan, D. and Henttonen, H. M. (1992). Marked Point Processes in Forest Statistics. \emph{Forest Science} 38(4): 806-824.\cr
	Penttinen, A. (2006). Statistics for Marked Point Patterns. in \emph{The Yearbook of the Finnish Statistical Society}. The Finnish Statistical Society, Helsinki: 70-91.
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}
\note{
	The function is computed using \code{\link{markcorrint}} in spatstat.
}

\seealso{
	\code{\link{Lmm.r}}, \code{\link{LmmEnvelope}}, \code{\link{markcorrint}}
}
\examples{

data(paracou16)
# Keep only 50% of points to run this example
X <- rthin(paracou16, 0.5)
plot(X)

# Calculate Kmm
r <- seq(0, 30, 2)
ActualValues.X <- Kmm.r(X, r)

# Plot
PlotResults(r, DivideByPiR2(ActualValues.X, r), ylab="Kmm / (pi R^2)", ReferenceValue=1)

}
