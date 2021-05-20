## Score equations for the marginalized two-part model

score = function(params, y, Z, X) {
  
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
  da = apply( as.vector(da)*Z, 2, sum)
  
  db = as.vector( { exp( -2 * s ) * xx + .5 } )
  db = apply(db * X[pos, 1:ncol(X), drop = FALSE], 2, sum)
  
  ds = sum( -1 + exp( -2 * s) * xxx^2  - xxx )

  -1*c( da, db, ds )
}

