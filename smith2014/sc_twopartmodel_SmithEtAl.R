# Example use of a two-part marginal model as in Smith (2014)

# libraries
library(stats4)

# data for testing
data(bioChemists, package = 'pscl')

# Initial parameters for the zero vs non-zero part: ---------------------------
logit_formula = I(art>0) ~ fem + mar + kid5 + phd + ment
logit = glm( as.formula(logit_formula), family=binomial(), data = bioChemists)
alpha_0 = coef(logit)

# Initial parameters for the continuous part: ---------------------------------
cont_formula = log(art) ~ fem + mar + kid5 + phd + ment
cont = lm( as.formula(cont_formula), data = bioChemists[bioChemists$art > 0, ])
beta_0 = coef(cont)
s0 = log( sigma(cont) )

# Design matrix for the zero-process: -----------------------------------------
Z = model.matrix(logit_formula, data = bioChemists)

# Design matrix for the continuous process: -----------------------------------
X = model.matrix(cont_formula, data = bioChemists)

# Set up a vector with all coefficients / parameters: -------------------------
params = vector(length = ncol(Z) + ncol(X) + 1, mode = 'numeric')
names(params) = c( paste0('zero_', colnames(Z)), 
                   paste0('mean_', colnames(X)), 
                   'log_sigma')
## Use initial values from above
params[ grep('zero_', names(params))] = alpha_0
params[ grep('mean_', names(params))] = beta_0
params["log_sigma"] = s0

# Log-liklihood for marginalize model from Smith (2014)
loglik = function(params, y, Z, X){
  
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
   sum( -1*log(1 + exp(p0)) ) +
   sum( p0[pos] - log( y[pos] ) - .5*log( 2*pi ) - s -
                  .5 / {exp(s)^2 } * { log( y[pos] ) + p0[pos] -
                      log( 1 + exp(p0[pos]) ) +
                      exp(s)^2*.5 - v[pos]  
                  }^2 
              )
}

# Test that loglik works properly
y = bioChemists$art 
loglik(params, y, Z, X)


# To pass this to stats4::mle, we need our function to only take the "params"
# argument.  
minus_loglik = function(theta, y, Z, X){
  -1*loglik(theta, y, Z, X)
}
## ? Why does this fail? 
#mle_fit = mle(minus_loglik, start = list(theta = params))

minus_loglik_grad = function(theta, y, Z, X){
  
  alpha = theta[grep('zero_', names(params))]
  beta = params[ grep('mean_', names(params))]
  sigma = exp(params['log_sigma'])

  Za = Z %*% alpha 
  Xb = X %*% beta
  
  pos = which(y > 0)
  term1 =  {1 - 
    1 / sigma^2 * { log(y) + Za - log( 1 + exp(Za)) + 1 / sigma^2 - Xb}
  }*{1 - plogis(Za)} 
  term2 = rep(0, length(y))
  term2[pos] = term1[pos]
  
  dl_dalpha = as.numeric( t(Z) %*% {-plogis( Za ) + term2} )

  term3 = {1 / sigma^2 *{ log(y) + Za  - log(1 + exp(Za)) - Xb} + .5}
  term3[is.infinite(term3)] = 0
  dl_dbeta = as.numeric( t(X) %*% term3 )
  
  term4 = {-log(sigma) -  
      .5 / sigma^2 * { log(y) + Za - log( 1 + exp(Za)) + sigma^2 - Xb}^2 - 
      { log(y) + Za - log( 1 + exp(Za)) + sigma^2 - Xb}
  }
  dl_dlogsigma = sum(term4[pos])

  -1 * c(dl_dalpha, dl_dbeta, dl_dlogsigma) 
  
}
length( minus_loglik_grad(params, y, Z, X) )

# We can maximize the liklihood manually: -------------------------------------
fit0 = optim(par = params, 
            f = minus_loglik,
            gr = minus_loglik_grad,
            method = 'BFGS',
            hessian = FALSE,  # We will compute after convergence
            y = y, Z = Z, X = X,
            control = list("maxit" = 1e4)
            )

fit1 = optim(par = params, 
             f = minus_loglik,
             hessian = FALSE,  # We will compute after convergence
             y = y, Z = Z, X = X,
             control = list("maxit" = 1e4)
)

fit0$value
fit1$value
cbind( fit0$par, fit1$par)
if( fit0$convergence != 0 ) {
  cat('Optimization may not have converged.\n')
}

fit = optim(par = fit0$par, 
             f = minus_loglik,
             hessian = TRUE,  # We will compute after convergence
             y = y, Z = Z, X = X,
             control = list("maxit" = 1e4)
)


# Estimated Coefficients
bhat = fit$par

# Estimated Variance-Covariance
v = solve( fit$hessian )

# Summary table 

# Standard error, z-score, and p-value.
std_error = sqrt(diag(v))
z = bhat / std_error
p = pmin( 2*{1 - pnorm(abs(z))}, 1)

# Stars denoting significance
stars = c('', '.', '*', '**', '***')[ sapply(p, function(x){
  1 + sum(x < c(.001, .01, .05, .1) )
  } ) ]
cbind("Estimate" = bhat, "Std. Error" = std_error, "z value" = z, "P(>|z|)" = p) 

# Compare to the hurdle model using a truncated Poisson for the count.
summary( pscl::hurdle(art ~ ., data = bioChemists) )

