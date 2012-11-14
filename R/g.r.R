g.r <-
function(X, r, ReferenceType="", NeighborType="") {
  # Eliminate erroneous configurations
  if ((ReferenceType=="" | NeighborType=="") & (ReferenceType!=NeighborType)) {
    stop("Either two or no point type must be specified.")
  }

  area <- area.owin(X$window)

  # g intra
  if (ReferenceType=="" & NeighborType=="") {
    Pairs <- closepairs(X, max(r))
    lambdaI <- lambdaJ <- X$n/area
  } else {
    # g intra for a single point type
    if (ReferenceType==NeighborType) {
      X.reduced <- X[X$marks$PointType==ReferenceType]
      Pairs <- closepairs(X.reduced, max(r))
      lambdaI <- lambdaJ <- X.reduced$n/area
    }
    # g inter
    if (ReferenceType!=NeighborType) {
      X.cross <- X[X$marks$PointType==ReferenceType]
      Y.cross <- X[X$marks$PointType==NeighborType]
      Pairs <- crosspairs(X.cross, Y.cross, max(r))
      lambdaI <- X.cross$n/area
      lambdaJ <- Y.cross$n/area
    }
  }

  # Adapted from pcf.ppp {spatstat}

  # Geometry
  XI <- ppp(Pairs$xi, Pairs$yi, window = X$window, check = FALSE)
  XJ <- ppp(Pairs$xj, Pairs$yj, window = X$window, check = FALSE)
  
  # Find the best breaks
  rmaxdefault <- rmax.rule("K", X$window, lambdaJ)
  breaks <- handle.r.b.args(window = X$window, rmaxdefault = rmaxdefault)
  rBest <- breaks$r
  denargs <- resolve.defaults(list(kernel = "epanechnikov",  n = length(rBest), from = 0, to = max(rBest)))

  # Edge-effect correction
  if (is.rectangle(X) |  is.polygonal(X)) {
    edgewt <- edge.Ripley(XI, matrix(Pairs$d, ncol = 1))
  } else {
    edgewt <- edge.Trans(XI, XJ, paired = TRUE)
  }

  # Estimate g  
  gEstimate <- sewpcf(Pairs$d, edgewt, denargs, lambdaI*lambdaJ*area)$g
  Approx <- approx(rBest, gEstimate, xout=r)$y
  # Eliminate inf values
  Approx[!is.finite(Approx)] <- NA
  
  # Return the values of g(r)
  return(Approx)
}
