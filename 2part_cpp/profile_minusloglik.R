# Function
minus_loglik = function(params, y, Z, X, w){
  
  # parameters
  alpha = params[ 1:ncol(Z) ] #grep('zero_', names(params))]
  beta = params[ {ncol(Z)+1}:{ncol(Z) + ncol(X)}]#grep('mean_', names(params))]
  
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
            exp(s)^2*.5 - v[pos]  }^2} 
    )
  
  return(-loglik)
}

# Test set up
source('~/github/2partmodels/smith2014/local/two_part_formuals.R')
data(bioChemists, package = 'pscl')
mf = two_part(art ~ phd + kid5 | fem + phd, data = bioChemists)
cont = lm(art ~ phd + kid5, data = bioChemists[bioChemists$art > 0,])
zero = glm( I(art > 0) ~ fem + phd, data = bioChemists, family = binomial())

X = mf$X; Y = mf$Y; Z = mf$Z
params = c(coef(zero), coef(cont), sigma(cont))
names(params) =  
  c( paste0('zero_', colnames(Z)), paste0('mean_', colnames(X)), 'log_sigma' )
w = rep(1, length(Y))

# profile
library(profvis)
w = rep(w, 100); 
Y = rep(Y, 100); dim(Y) = c(length(Y), 1)
X = rep(X, 100); dim(X) = c(length(Y), length(coef(cont)))
Z = rep(Z, 100); dim(Z) = c(length(Y), length(coef(zero)))

microbenchmark::microbenchmark( minus_loglik(params, Y, Z, X, w), times = 10 )
profvis({
  minus_loglik(params, Y, Z, X, w)
})
