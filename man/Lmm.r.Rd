\name{Lmm.r}
\alias{Lmm.r}
\title{
  Estimation of the Lmm function
}
\description{
  Estimates the Lmm function
}
\usage{
Lmm.r(X, r, ReferenceType = "")
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
	Lmm is the normalized version of Kmm: Lmm = sqrt(Kmm/pi)-r
}
\value{
	A vector containing L(r) values.
}
\references{
	Penttinen, A., Stoyan, D. and Henttonen, H. M. (1992). Marked Point Processes in Forest Statistics. Forest Science 38(4): 806-824.\cr
	Espa, G., Giuliani, D. and Arbia, G. (2010). Weighting Ripley's K-function to account for the firm dimension in the analysis of spatial concentration. Discussion Papers, 12/2010. Universita di Trento, Trento: 26.	
}
\author{
	Eric Marcon <Eric.Marcon@ecofog.gf>
}


\seealso{
	\code{\link{Kmm.r}}, \code{\link{LmmEnvelope}}
}
\examples{

data(paracou16)
# Keep only 50% of points to run this example
X <- rthin(paracou16, 0.5)
plot(X)

# Calculate Lmm
r <- seq(0, 30, 2)
ActualValues.X <- Lmm.r(X, r)

# Plot
PlotResults(r, ActualValues.X, ylab="Lmm", ReferenceValue=0)

}
