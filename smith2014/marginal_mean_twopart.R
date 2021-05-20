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

marginal_mean_twopart = function(formula, data, 
                                 weights = NULL,
                                 optim.control = NULL, 
                                 method = 'BFGS'){
  # Inputs: 
  #   formula - a formula of the form response ~ z1 + z2 ... | x1 + x2 ...
  #             where z1/x1 are possibly the same.  A model with no bar "|" will
  #             use the same terms for both the logit and log-normal models. 
  #   data - a data frame where the terms in the formula can be found
  #   weights - an optional vector of weights for each case
  
  ## First, split the formula and form the needed model matrices: ------------- 
  ## extracted from pscl::zero_inflation, identifies "|" in formulas. 
  cl <- match.call() 
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)#;  return(mf)
  
  m <- match(c("formula", "data", "subset", "na.action", "weights", 
               "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  #! Recently added, check here if this breaks. 
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")
  
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
  # w - is a case weight or the constant 1, if NULL

  w = weights
  if ( is.null(w) ) w = rep(1, nrow(X) )
  
  # Get initial starting values using distinct models
  # Initial parameters for the zero vs non-zero part: -------------------------
  Y1 = 1*{Y > 0}
  
  if ( length(w) == nrow(Z) ) {
    suppressWarnings({
      logit = glm.fit(x = Z, y = Y1, family=binomial(), weights = w)
    })
  }
  else {
    logit = glm.fit(x = Z, y = Y1, family=binomial())
  }
    
  alpha_0 = coef(logit)
  
  # Initial parameters for the continuous part: ------------------------------
  #! Use log Y in the initial fit to start closer to solution.
  if ( length(w) == sum( Y1 == 1 ) ) {
    cont = lm.wfit(X[ Y1 == 1, 1:ncol(X) , drop = FALSE], log( Y[Y1 == 1] ),
                   w = w)    
  } else {
    cont = lm.fit(X[ Y1 == 1, 1:ncol(X) , drop = FALSE], log( Y[Y1 == 1] ) )
  }

  beta_0 = coef(cont)
  
  r =  log( Y[ Y1 == 1 ] ) - X[ Y1 == 1, 1:ncol(X), drop = FALSE] %*% beta_0
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
  minus_loglik = function(params, y, Z, X, w){
    
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
    loglik = sum( -w * log( 1 + exp(p0) ) ) +
      sum( w[pos]*{p0[pos] - log( y[pos] ) - .5*log( 2*pi ) - s -
             .5 / {exp(s)^2 } * { log( y[pos] ) + p0[pos] -
                 log( 1 + exp(p0[pos]) ) +
                 exp(s)^2*.5 - v[pos]  
             }^2} 
       )
    
    return(-loglik)
  }
  #minus_loglik(params, y, Z, X, w)  
  
  ## Score equations for the marginalized two-part likelihoo above
  score = function(params, y, Z, X, w) {
    
    # sZ = apply(Z, 2, sum)
    # sX = apply(X, 2, sum)
    
    # parameters
    alpha = params[ grep('zero_', names(params))]
    beta = params[ grep('mean_', names(params))]
    
    # compute matrix products
    p0 = Z %*% alpha
    v = X %*% beta
    
    pos = which(y > 0)
    n1 = sum(pos)
    s = params["log_sigma"]
    
    xx = log(y[pos]) + p0[pos] - log( 1 + exp(p0[pos]) ) - v[pos]
    xxx = xx + .5 * exp( 2 * s)
    
    da =  -plogis( p0 )
    da[pos] = da[pos] + { 1 - exp( -2 * s ) * xxx * plogis( -p0[pos] ) }
    da = apply( w * as.vector(da) * Z, 2, sum)
    
    db = w[pos]*as.vector( { exp( -2 * s ) * xx + .5 } )
    db = apply( db * X[pos, 1:ncol(X), drop = FALSE], 2, sum)
    
    ds = sum( w[pos] * {-1 + exp( -2 * s) * xxx^2  - xxx} )
    
    -1*c( da, db, ds )
  }
  #score(params, y, Z, X, w = runif(length(y)))
  
 
  # Directly optimize: --------------------------------------------------------
  fit1 = optim(par = params, 
               f = minus_loglik,
               gr = score,
               method = method,
               hessian = TRUE,  
               y = Y, Z = Z, X = X, w = w, 
               control = optim.control
  )
  
  #! Also do some NA checking up front.
  #return( fit1 )

  if ( fit1$convergence != 0 )  cat(fit1a$message,'\n')

  if ( nrow(X) != nrow(Z) ) warning("nrow(X) and nrow(Z) are not equal.\n")
  
  out = list( coefficients = fit1$par, 
              logLik = fit1$value,
              convergence = fit1$convergence,
              termsX = mtX,
              termsZ = mtZ,
              call = cl,
              vcov = solve( fit1$hessian ),
              n = min( nrow(X), nrow(Z) ),
	      logit = logit,
	      cont = cont
              )
  class(out) = c('mm2p')

  out
}

# Some generic methods: -------------------------------------------------------
print.mm2p = function( obj ) {
  getS3method("print", "lm")(obj)
  cat('n:\n', obj$n, '\n')
}

summary.mm2p = function( obj ) {
  
  # Summary table 
  
  # Standard error, z-score, and p-value.
  v = obj$vcov
  std_error = sqrt(diag(v))
  z = obj$coefficients / std_error
  p = pmin( 2*{1 - pnorm(abs(z))}, 1)
  
  # Stars denoting significance
  stars = c('', '.', '*', '**', '***')[ sapply(p, function(x){
    1 + sum(x < c(.001, .01, .05, .1) )
  } ) ]
  
  coefficients = 
    cbind( "Estimate" = obj$coefficients, 
           "Std. Error" = std_error, 
           "z value" = z, 
           "P(>|z|)" = p
    )
  
  # 
  out = list( call = obj$call, coefficients = coefficients)
  class( out ) = c('summary_mm2p')
  out
}

print.summary_mm2p = function( obj ) {
  # Stars denoting significance
  
  p = obj$coefficients[,"P(>|z|)"] 
  stars = c('', '.', '*', '**', '***')[ sapply(p, function(x){
    1 + sum(x < c(.001, .01, .05, .1) )
  } ) ]
  
  out = as.data.frame(obj$coefficients)
  out$' ' = stars
  
  print(out)
}

predict.mm2p = 
 function(obj, newdata, 
          type = c('marginal_mean', 'conditional_mean', 'prob'),
          se.method = c('none', 'dm') ) { 
   # marginal_mean: predicts the marginal mean for each row of newdata
   # condtional_mean: predicts the mean conditional on being non-zero
   # prob: predicts the probability that the case has a non-zero respons
   # se.method: the method used to compute standard errors, the default "none"
   # results in no standard errors being returned
            
  type = match.arg(type, c('marginal_mean', 'conditional_mean', 'prob') )
  se.method = match.arg(se.method, c('none', 'dm'))
  
  if ( !{ type %in% c('marginal_mean', 'conditional_mean', 'prob')} ) {
    stop("please specify a valid type")
  }
  if ( type == 'marginal_mean' ) {
     X = model.matrix( update( obj$termsX, NULL ~ .), data = newdata )
     yhat = exp( 
       X %*% obj$coefficients[ grep('^mean_', names(obj$coefficients) ) ] 
     )
     
     # delta method standard errors
     if ( se.method == 'dm' ) {
       mmind = grep('^mean_', colnames( obj$vcov ) )
       v = obj$vcov[mmind, mmind]
       yhat_se = sqrt( diag( X %*% v %*% t(X) ) ) 
     }

  }
  
  if ( type == 'prob' ) {
    Z = model.matrix( update( obj$termsZ, NULL ~ .), data = newdata )
    yhat = plogis( 
      Z %*% obj$coefficients[ grep('^zero_', names(obj$coefficients) ) ] 
    )
  }
  
  if ( type == 'conditional_mean' ) {
    Z = model.matrix( update( obj$termsZ, NULL ~ .), data = newdata )
    phat = plogis( 
      Z %*% obj$coefficients[ grep('^zero_', names(obj$coefficients) ) ] 
    )
    
    X = model.matrix( update( obj$termsX, NULL ~ .), data = newdata )
    yhat = exp( 
      X %*% obj$coefficients[ grep('^mean_', names(obj$coefficients) ) ] 
    )
    
    yhat = yhat / phat
  }
  
  if ( se.method == 'none' ){
    return( as.vector( yhat ) )
  } else {
    return( data.frame( yhat = yhat, se = yhat_se) )
  }

}

# Gradient/score function for the likelihood above

## For testing internally.
if ( FALSE ) {
  source('~/github/2partmodels/smith2014/local/two_part_formuals.R')
  data(bioChemists, package = 'pscl')
  mf = two_part(art ~ phd + kid5 | fem + phd, data = bioChemists)
  cont = lm(art ~ phd + kid5, data = bioChemists[bioChemists$art > 0,])
  zero = glm( I(art > 0) ~ fem + phd, data = bioChemists, family = binomial())
  sigma(cont)
  X = mf$X; Y = mf$Y; Z = mf$Z
  params = c(coef(zero), coef(cont), sigma(cont))
  names(params) =  c( paste0('zero_', colnames(Z)), 
                      paste0('mean_', colnames(X)), 
                      'log_sigma')
  
  #score(params, Y, Z, X)
  #params = params; y = Y; 
  ## Example of how to run
  #fit1 = marginal_mean_twopart( art ~ phd + kid5 | fem + phd, data = bioChemists,
  #                              optim.control = list(maxit = 1e5), 
  #                              method = 'Nelder-Mead')
  #fit1$convergence
  #fit1$counts
  
  #fit2 = marginal_mean_twopart( art ~ phd + kid5 | fem + phd, data = bioChemists,
  #                              optim.control = list(maxit = 1e5),
  #                              method = 'BFGS')
  #fit2$convergence
  #fit2$counts
  #cbind(fit1$par, fit2$par)
}
