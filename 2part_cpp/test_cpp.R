library(Rcpp)

sourceCpp("./mm2p_loglik.cpp")
sourceCpp("./log_lik_sum.cpp")

.5 * log(2 * pi) + s
v3 = log( y[pos] ) + p0[pos] - log( 1 + exp(p0[pos]) ) +
  .5 * exp(2 * s) - #v2
  v[pos] 
