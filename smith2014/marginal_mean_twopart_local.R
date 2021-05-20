# Function to fit a marginalized two-part (hurdle) model a la Smith et al 2014.
#
# The model consists of a logistic regression model for zero vs non-zero costs
# and a log Normal specification for the magnitude of the non-zero costs.
#
# This is a separable two-part model, and could, in theory, be fit as two 
# distinct models. However, rather than marginalizing these disinct models, 
# the joint-likelihood is parameterized in terms of the marginalized mean and
# maximized.
#
# The formula specification uses "|" on the right-hand side to allow different
# covariates to be used for the continuous and logistic models.  For details,
# see the help('zeroinfl', package = 'pscl')
#
# Author: James Henderson, PhD (jbhender@umich.edu)
# Date: June 28, 2019

marginal_mean_twopart = function(formula, data, optim.control = NULL){
  
  ## First, split the formula and form the needed model matrices: ------------- 
  ## extracted from pscl::zero_inflation, identifies "|" in formulas. 
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)#;  return(mf)
  
  m <- match(c("formula", "data", "subset", "na.action", "weights", 
               "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  if ( length(formula[[3]]) > 1 && 
       identical(formula[[3]][[1]], as.name("|")) ) {
    ff <- formula
    formula[[3]][1] <- call("+")
    mf$formula <- formula
    ffc <- . ~ .
    ffz <- ~.
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
  } else {
    ffz <- ffc <- ff <- formula
    ffz[[2]] <- NULL
  }
  if (inherits(try(terms(ffz), silent = TRUE), "try-error")) {
    ffz <- eval(parse(text = sprintf(paste("%s -", deparse(ffc[[2]])), 
                                     deparse(ffz))))
  }
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  mtX <- terms(ffc, data = data)
  X <- model.matrix(mtX, mf)
  mtZ <- terms(ffz, data = data)
  mtZ <- terms(update(mtZ, ~.), data = data)
  Z <- model.matrix(mtZ, mf)
  Y <- model.response(mf, "numeric")
  ## end code taken from pscl::zero_inflation
  
  # Notes: --------------------------------------------------------------------
  # Y - is the response (i.e. costs) 
  # X - is the model matrix for the continuous log-normal model
  # Z - is the model matrix for the logistic model
  
  # Get initial starting values using distinct models
  # Initial parameters for the zero vs non-zero part: -------------------------
  Y1 = 1*{Y > 0}
  pos = which(Y1 > 0)
  logit = glm.fit(x = Z, y = Y1, family=binomial())
  alpha_0 = coef(logit)
  
  # Initial parameters for the continuous part: ------------------------------
  cont = lm.fit(X[ pos, ], Y[pos]) 
  beta_0 = coef(cont)
  r = Y[ pos ] - X[ pos, ] %*% beta_0
  s0 = log( sqrt( sum(r^2) / {length(r) - ncol(X)}) )

  # Set up a vector with all coefficients / parameters: -----------------------
  params = vector(length = ncol(Z) + ncol(X) + 1, mode = 'numeric')
  names(params) = c( paste0('zero_', colnames(Z)), 
                     paste0('mean_', colnames(X)), 
                     'log_sigma')
  ## Use initial values from above
  params[ grep('zero_', names(params))] = alpha_0
  params[ grep('mean_', names(params))] = beta_0
  params["log_sigma"] = s0

  # Log-liklihood for marginalize model from Smith (2014): --------------------
  minus_loglik = function(params, y, Z, X){
    
    # parameters
    alpha = params[ grep('zero_', names(params))]
    beta = params[ grep('mean_', names(params))]
    
    # compute matrix products
    p0 = Z %*% alpha
    v = X %*% beta
    
    # non-zero locations in y
    pos = which(y > 0)
    n1 = length(pos)
    s = params["log_sigma"]
    
    # log likelihood
    loglik = 
      sum( -1*log(1 + exp(p0)) ) +
      sum( p0[pos] - log( y[pos] ) - .5*log( 2*pi ) - 
           s - .5 / {exp(s)^2 } * 
           { log( y[pos] ) + p0[pos] - log( 1 + exp(p0[pos]) ) + 
             exp(s)^2*.5 - v[pos]  
           }^2 
       )
    
    return(-loglik)
  }
 
  # Directly optimize: --------------------------------------------------------
  fit1 = optim(par = params, 
               f = minus_loglik,
               hessian = TRUE,  
               y = Y, Z = Z, X = X,
               control = optim.control
  )
  
  return( fit1 )
}

## For testing internally.
# 
# mf = two_part(art ~ phd + kid5 | fem + phd, data = bioChemists)
# cont = lm(art ~ phd + kid5, data = bioChemists[bioChemists$art > 0,])
#sigma(cont)
# X = mf$X; Y = mf$Y; Z = mf$Z

## Example of how to run
fit1a = marginal_mean_twopart( art ~ phd + kid5 | fem + phd, data = bioChemists,
                              optim.control = list(maxit = 1e5))
fit1a$convergence
names(fit1a)
fit1a$counts

#fit1a$par == fit1b$par
