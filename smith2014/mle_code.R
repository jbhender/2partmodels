# our values
minuslogl = minus_loglik
start = list(theta = params) 
method = "BFGS"
fixed = list()

# function "mle"
n <- names(fixed)
fullcoef <- formals(minuslogl)
if (any(!n %in% names(fullcoef))) 
  stop("some named arguments in 'fixed' are not arguments to the supplied log-likelihood")
fullcoef[n] <- fixed
if (!missing(start) && (!is.list(start) || is.null(names(start)))) 
  stop("'start' must be a named list")
start[n] <- NULL
# ??
#start <- sapply(start, eval.parent)
nm <- names(start)
oo <- match(nm, names(fullcoef))
if (anyNA(oo)) 
  stop("some named arguments in 'start' are not arguments to the supplied log-likelihood")
start <- start[order(oo)]
nm <- names(start)
f <- function(p) {
  l <- as.list(p)
  names(l) <- nm
  l[n] <- fixed
  do.call("minuslogl", l)
}
# Fails during these lines
oout <- if (length(start)) 
  optim(start, f = f, method = method, hessian = TRUE)
  else list(par = numeric(), value = f(start))
coef <- oout$par
vcov <- if (length(coef)) 
  solve(oout$hessian)
else matrix(numeric(), 0L, 0L)
min <- oout$value
fullcoef[nm] <- coef
new("mle", call = call, coef = coef, fullcoef = unlist(fullcoef), 
    vcov = vcov, min = min, details = oout, minuslogl = minuslogl, 
    nobs = if (missing(nobs)) 
      NA_integer_
    else nobs, method = method)
