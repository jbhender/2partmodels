# Illustrations of the log-normal distribution


n = 100
z = rnorm(n, 0, 1)
mu = 1; s = 1


x = exp(mu + s*z)

muhat = mean( log(x) )
shat = sd( log(x) )*sqrt ( {n - 1} / n )

m = mean( x )
v = var(x)

log( m / sqrt(1 + v / m^2 ) ); muhat
exp( muhat + shat^2 / 2); m
( exp( shat^2 ) - 1 ) * exp( 2*muhat + shat^2 ); v
log( 1 + v / m^2); shat

library(EnvStats)
elnormAlt(x)
