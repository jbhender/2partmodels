# Trouble-shoot Erin Ice's Example for Smith et al's 2014 marginalized two-part
# model code.
#
# Date: July 18, 2019
# Author: James Henderson, PhD (jbhender@umich.edu)

# libraries: ------------------------------------------------------------------
library(tidyverse); library(foreign); library(EnvStats)

# source code: ----------------------------------------------------------------
source('./marginal_mean_twopart.R')

# data for testing: -----------------------------------------------------------
atus_long <- read.dta("~/Dropbox/marginal_mean_twopart_ice/atus_longmerge_r.dta")
atus_long <- atus_long %>% filter(year >= 2005)

# Fit the model as in the original function, using Nelder-Mead to optimize: ---
if ( FALSE ) {
 tm0 = system.ti1me({
 fit1a = marginal_mean_twopart(
   totalcare ~ factor(year) + male*ba + age + raceeth + weekend | 
              factor(year) + male*ba + age + raceeth + weekend, 
              data = atus_long, method = 'Nelder-Mead',
              optim.control = list(maxit = 1e5))
 })
}
#fit1a$counts
if( fit1a$convergence != 0 ) cat('Model did not converge.\n')

# Optimize using the BFGS algorithm and the (corrected) gradient: -------------
tm1 = system.time({
  fit1b = marginal_mean_twopart(
  totalcare ~ factor(year) + male*ba + age + raceeth + weekend | 
    factor(year) + male*ba + age + raceeth + weekend, 
  data = atus_long, method = 'BFGS',
  optim.control = list(maxit = 1e5))
})


# You can compare tm0 and tm1 to see if one is quicker.  
# Estimated Coefficients
bhat = fit1b$par

# Estimated Variance-Covariance
v = solve( fit1b$hessian )

# Summary table 

# Standard error, z-score, and p-value.
std_error = sqrt(diag(v))
z = bhat / std_error
p = pmin( 2*{1 - pnorm(abs(z))}, 1)


# Stars denoting significance
stars = c('', '.', '*', '**', '***')[ sapply(p, function(x){
  1 + sum(x < c(.001, .01, .05, .1) )
} ) ]

data.frame(
  "Estimate" = bhat, "Std. Error" = std_error, "z value" = z, "P(>|z|)" = p, 
  "Stars" = stars )
 

table1 <- data.frame(
  "Estimate" = round(bhat, 3), 
  "Std. Error" = round(std_error, 2), 
  " " = stars )


# Non-marginal mean, conditional on totalcare > 0: ----------------------------
fit_pos = lm( log( totalcare) ~ 
                factor(year) + male*ba + age + raceeth + weekend,
              data = atus_long %>% filter(totalcare > 0) )
summary(fit_pos)

fit_zero = glm( I( totalcare > 0 ) ~ 
                factor(year) + male*ba + age + raceeth + weekend,
              data = atus_long, family = binomial( link = 'logit') )

summary(fit_zero)

atus_long %>%
  group_by(weekend) %>%
  summarize( phat = mean( totalcare > 0 ), 
             xbar = mean( totalcare[totalcare > 0] ) 
  )
  
table1

## Table 1
fit0a = marginal_mean_twopart(
  totalcare ~ factor(year) | factor(year), 
  data = atus_long, method = 'BFGS',
  optim.control = list(maxit = 1e3))

fit0b = marginal_mean_twopart(
  totalcare ~ factor(year) | factor(year), 
  data = atus_long, method = 'Nelder-Mead',
  optim.control = list(maxit = 1e3))

round( cbind(fit0a$par, fit0b$par), 3 )

## Intercept only model
fit2a = marginal_mean_twopart(
  totalcare ~ 1 | 1, 
  data = atus_long, method = 'BFGS',
  optim.control = list(maxit = 1e3))

fit2b = marginal_mean_twopart(
  totalcare ~ 1 | 1, 
  data = atus_long, method = 'Nelder-Mead',
  optim.control = list(maxit = 1e3))

round( cbind(fit1a$par, fit1b$par), 3 )
exp( fit1a$par[2] ) / plogis( fit1a$par[1] )
exp( fit1b$par[2] ) / plogis( fit1b$par[1] )

bhat = fit0$par

# Estimated Variance-Covariance
v = solve( fit0$hessian )

# Summary table 

# Standard error, z-score, and p-value.
std_error = sqrt(diag(v))
z = bhat / std_error

data.frame( estimate = round( bhat, 2), ncol = 1)

atus_long %>% group_by(year) %>% summarize( m = mean(totalcare) )
atus_long %>% summarize( z = mean(totalcare > 0),
                         m = mean(totalcare), 
                         mpos = mean( totalcare[ totalcare > 0]),
                         mlpos = mean( log( totalcare[ totalcare > 0] ) ),
                         sigma = sd( log( totalcare[ totalcare > 0])),
                         v = var( totalcare[totalcare > 0])
                         )

fit2z = glm( I(totalcare > 0) ~ 1, data = atus_long, family = binomial())
summary(fit2z)

fit2p = lm( log(totalcare) ~ 1, data = atus_long %>% filter(totalcare > 0))
summary(fit2p)

n = sum( atus_long$totalcare > 0 )
s = sigma(fit2p)*{ n - 1} / n
plogis( coef(fit2z) ) *  exp( coef(fit2p) + s^2/2 ) 

log( 1 + 16134.05 / {114.4857^2} )
(exp( s^2) - 1)*exp(2*coef(fit2p) + s^2)

elnormAlt( atus_long$totalcare[ atus_long$totalcare > 0])
.46 * 149

# 

m1 = glm( totalcare ~ factor(year) + male*ba + age + raceeth + weekend,
          family = gaussian( link = 'log' ),
     data = atus_long %>% filter(totalcare > 0) )

m2 = glm( totalcare ~ factor(year) + male*ba + age + raceeth + weekend,
          family = Gamma( link = 'log' ),
          data = atus_long %>% filter(totalcare > 0) )

m3 = glm( totalcare ~ factor(year) + male*ba + age + raceeth + weekend,
          family = Gamma( link = 'inverse' ),
          data = atus_long %>% filter(totalcare > 0) )

m4 = glm( totalcare ~ factor(year) + male*ba + age + raceeth + weekend,
          family = Gamma( link = 'identity' ),
          data = atus_long %>% filter(totalcare > 0) )

BIC(fit_pos, m1, m2, m3, m4)

atus_long_nonzero = atus_long %>%
  filter(totalcare > 0) %>%
  mutate( mu1 = predict(m1, type = 'response'), 
          mu2 = predict(m2, type = 'response'),
          mu3 = predict(m3, type = 'response'),
          mu4 = predict(m4, type = 'response')
          )

atus_long_nonzero %>%
  summarize( m1 = mean(mu1), m2 = mean(mu2), m3 = mean(mu3), m4 = mean(mu4))
