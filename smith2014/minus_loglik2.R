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

minus_loglik2 = function(params, ypos, Z, X, ypos_logsum = 0){
  ## ypos - values of non-zero values
  ##    Z - model matrix for all terms, zero and nonzero
  ##    X - model matrix for positive terms only
  ## ypos_logsum - we can omit this and add it back in after optimization
  
  # parameters 
  ## Index rather than name, skip intermediate assignment
  
  # compute matrix products
  p0 = Z %*%  params[ 1:ncol(Z) ]
  v = X %*% params[  { ncol(Z) + 1 }:{ ncol(Z) + ncol(X) }]
  s = params["log_sigma"]
  p0pos = p0[ {nrow(Z) - nrow(X) + 1}:nrow(Z) ]
    
  # non-zero locations in y
  ## reorder to be in blocks
  #pos = which(y > 0)
  n1 = nrow(X)
  
  
  # log likelihood
  ## Don't depend on parameters, ignore until after optimization.
  #- log( y[pos] ) - .5*log( 2*pi ) - 
  
  ## refactor to be more efficient
  -1*{
   sum( -1*log(1 + exp(p0)) ) +
   sum(p0pos) -
   n1 * s + n1 * .125 * exp(2*s)  -  
   1 / {2*exp(s)^2 } * 
    sum( { log( ypos ) +
        p0pos - log( 1 + exp(p0pos) ) - v }^2 ) 
  }
}

## Test that the two functions are the same, up to the additive constant:
data(bioChemists, package = 'pscl')
source('~/github/2partmodels/smith2014/local/two_part_formuals.R')
mf = two_part(art ~ phd + kid5 | fem + phd, data = bioChemists)
cont = lm(art ~ phd + kid5, data = bioChemists[bioChemists$art > 0,])
zero = glm( I(art>0) ~ fem + phd, data = bioChemists, family = binomial())

#sigma(cont)
X = mf$X; Y = mf$Y; Z = mf$Z

#params = fit1a$par
minus_loglik(params, Y, Z, X)

ypos = Y[ Y > 0]
Z2 = Z[ c(which(Y == 0), which(Y != 0)), ]
cnst = {1 + .5*exp(params['log_sigma']*2)}*sum(log(ypos)) + length(ypos)*.5*log(2*pi)

minus_loglik2(params, ypos, Z2, X = X[Y > 0, ]) + cnst
