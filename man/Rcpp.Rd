\name{CountNbd}
\alias{CountNbd}
\alias{CountNbdCC}
\alias{DistKd}
\alias{CountNbdKd}
\title{
  Auxiliary functions to count point neighbors
}
\description{
  C++ routines used for fast count of neighbors.
}
\usage{
CountNbd(Rr, Rx, Ry, RWeight, RNbd, RIsReferenceType, RIsNeighborType)
CountNbdCC(Rr, Rx, Ry, RWeight, RNbd, RIsReferenceType, RIsNeighborType)
DistKd(Rx, Ry, RPointWeight, RWeights, RDist, RIsReferenceType, RIsNeighborType)
CountNbdKd(Rr, Rx, Ry, RWeight, RNbd, RIsReferenceType, RIsNeighborType)
}
\arguments{
  \item{Rr}{
    A vector of distances.
  }
  \item{Rx}{
    A vector of \emph{x} coordinates.
  }
  \item{y}{
    A vector of \emph{y} coordinates.
  }
  \item{Weight}{
    A vector of point weights.
  }
  \item{PointWeight}{
    A vector of point weights.
  }
  \item{Nbd}{
    A matrix to store the weight of neighbors.
  }
  \item{Weights}{
    A vector of point-pair weights.
  }
  \item{Dist}{
    A vector of point-pair distances.
  }
  \item{IsReferenceType}{
    A boolean vector to recognize reference points.
  }
  \item{IsNeighborType}{
    A boolean vector to recognize neighbor points.
  }
}
\details{
  These routines are called internally.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\keyword{internal}