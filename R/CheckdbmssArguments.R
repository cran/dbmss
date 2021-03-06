CheckdbmssArguments <-
function() {

  # Verify that the package is attached
  if (! "dbmss" %in% .packages()) {
    warning("Function arguments cannot be checked because the dbmss package is not attached. Add CheckArguments=FALSE to suppress this warning or run library('dbmss')")
    return (TRUE)
  }
  # Get the list of arguments of the parent function
  ParentFunction <- sys.call(-1)[[1]]
  # If apply() or similar was used, the function name is not in ParentFunction: sys.call(-1)[[1]] returns "FUN"
  if (ParentFunction == "FUN") {
    warning("Function arguments cannot be checked, probably because you used apply(). Add CheckArguments=FALSE to suppress this warning.")
    return (TRUE)
  }
  
  # Find the arguments. match.fun does not work with dbmss::function
  # as.character creates a vector. The name of the function is the last item
  ParentFunction_split <- as.character(ParentFunction)
  ParentFunctionNoNS <- ParentFunction_split[length(ParentFunction_split)]
  Args <- formals(match.fun(ParentFunctionNoNS))
  
  ErrorFunction <- paste("Error in ", ParentFunctionNoNS, ":")
  
  ErrorMessage <- function(Message, Argument) {
    cat(deparse(substitute(Argument)), "cannot be:\n")
    print(utils::head(Argument))
    cat(paste(ErrorFunction, Message, "\n"))
    stop("Check the function arguments.", call. = FALSE)
  }
  

  # Get the point pattern or the Dtable
  X <- eval(expression(X), parent.frame())
  # X 
  if (!is.na(names(Args["X"]))) {
    if (!(inherits(X, "wmppp") | (inherits(X, "Dtable"))))
      ErrorMessage("X must be of class wmppp or Dtable", X)
  }

  # r
  if (!is.na(names(Args["r"]))) {
    r <- eval(expression(r), parent.frame())
    if (!is.null(r)) {
      if (!is.numeric(r) && !is.vector(r)) 
        ErrorMessage("r must be a numeric vector", r)
      if (length(r) < 2) 
        ErrorMessage(paste("r has length", length(r), "- must be at least 2"), r)
      if (r[1] != 0) 
        ErrorMessage("First r value must be 0", r)
      if (any(diff(r) <= 0)) 
        ErrorMessage("successive values of r must be increasing", r)  
    }
  }
    
  # ReferenceType 
  if (!is.na(names(Args["ReferenceType"]))) {
    ReferenceType <- eval(expression(ReferenceType), parent.frame())
    if (ReferenceType!="" & !ReferenceType %in% X$marks$PointType)
      ErrorMessage("ReferenceType must be a point type of the point pattern", ReferenceType)
  }
  # NeighborType 
  if (!is.na(names(Args["NeighborType"]))) {
    NeighborType <- eval(expression(NeighborType), parent.frame())
    if (NeighborType!="" & !NeighborType %in% X$marks$PointType)
      ErrorMessage("NeighborType must be a point type of the point pattern", NeighborType)
  }
  # Cases 
  if (!is.na(names(Args["Cases"]))) {
    Cases <- eval(expression(Cases), parent.frame())
    if (!Cases %in% X$marks$PointType)
      ErrorMessage("Cases must be a point type of the point pattern", Cases)
  }
  # Controls 
  if (!is.na(names(Args["Controls"]))) {
    Controls <- eval(expression(Controls), parent.frame())
    if (!is.null(Controls)) {
      if (!(Controls %in% X$marks$PointType))
        ErrorMessage("Controls must be a point type of the point pattern", Controls)
    }
  }
  
  # CaseControl 
  if (!is.na(names(Args["CaseControl"]))) {
    CaseControl <- eval(expression(CaseControl), parent.frame())
    if (!is.logical(CaseControl))
      ErrorMessage("CaseControl must be TRUE or FALSE", CaseControl)
  }
  # Intertype 
  if (!is.na(names(Args["Intertype"]))) {
    Intertype <- eval(expression(Intertype), parent.frame())
    if (!is.logical(Intertype))
      ErrorMessage("Intertype must be TRUE or FALSE", Intertype)
  }
  # Weighted 
  if (!is.na(names(Args["Weighted"]))) {
    Weighted <- eval(expression(Weighted), parent.frame())
    if (!is.logical(Weighted))
      ErrorMessage("Weighted must be TRUE or FALSE", Weighted)
  }
  # Original 
  if (!is.na(names(Args["Original"]))) {
    Original <- eval(expression(Original), parent.frame())
    if (!is.logical(Original))
      ErrorMessage("Original must be TRUE or FALSE", Original)
  }

  # lambda
  if (!is.na(names(Args["lambda"]))) {
    lambda <- eval(expression(lambda), parent.frame())
    if (!is.null(lambda)) {
      if (!inherits(lambda, "im") & !is.numeric(lambda))
        ErrorMessage("lambda must be an image of class im or a numeric vector", lambda)
    }
  }
  
  # NumberOfSimulations 
  if (!is.na(names(Args["NumberOfSimulations"]))) {
    NumberOfSimulations <- eval(expression(NumberOfSimulations), parent.frame())
    if (!is.numeric(NumberOfSimulations))
      ErrorMessage("NumberOfSimulations must be a number", NumberOfSimulations)
    if (NumberOfSimulations <= 0)
      ErrorMessage("NumberOfSimulations must be positive", NumberOfSimulations)
  }
  
  # Alpha 
  if (!is.na(names(Args["Alpha"]))) {
    Alpha <- eval(expression(Alpha), parent.frame())
    if (!is.numeric(Alpha))
      ErrorMessage("Alpha must be a number", Alpha)
    if (Alpha < 0)
      ErrorMessage("Alpha must be positive", Alpha)   
  }
  
  # alpha 
  if (!is.na(names(Args["alpha"]))) {
    alpha <- eval(expression(alpha), parent.frame())
    if (!is.numeric(alpha))
      ErrorMessage("alpha must be a number", alpha)
    if (alpha < 0)
      ErrorMessage("alpha must be positive", alpha)
    if (alpha > 1)
      ErrorMessage("alpha must be less than or equal to 1", alpha)
  }
  
  # Adjust 
  if (!is.na(names(Args["Adjust"]))) {
    Adjust <- eval(expression(Adjust), parent.frame())
    if (!is.numeric(Adjust))
      ErrorMessage("Adjust must be a number", Adjust)
    if (Adjust<=0)
      ErrorMessage("Adjust must be strictly positive", Adjust)
  }
  
  # Approximate 
  if (!is.na(names(Args["Approximate"]))) {
    Approximate <- eval(expression(Approximate), parent.frame())
    if (!is.numeric(Approximate))
      ErrorMessage("Approximate must be a number", Approximate)
    if (Approximate < 0)
      ErrorMessage("Approximate must be positive", Approximate)
  }

  # StartFromMinR 
  if (!is.na(names(Args["StartFromMinR"]))) {
    StartFromMinR <- eval(expression(StartFromMinR), parent.frame())
    if (!is.logical(StartFromMinR))
      ErrorMessage("StartFromMinR must be TRUE or FALSE", StartFromMinR)
  }
  
  # Individual 
  if (!is.na(names(Args["Individual"]))) {
    Individual <- eval(expression(Individual), parent.frame())
    if (!is.logical(Individual))
      ErrorMessage("Individual must be TRUE or FALSE", Individual)
  }
  
  # Precision 
  if (!is.na(names(Args["Precision"]))) {
    Precision <- eval(expression(Precision), parent.frame())
    if (!is.numeric(Precision))
      ErrorMessage("Precision must be a number", Precision)
    if (Precision < 0)
      ErrorMessage("Precision must be positive", Precision)
  }

  # show.window 
  if (!is.na(names(Args["show.window"]))) {
    show.window <- eval(expression(show.window), parent.frame())
    if (!is.logical(show.window))
      ErrorMessage("show.window must be TRUE or FALSE", show.window)
  }
  
  return (TRUE)
}
    