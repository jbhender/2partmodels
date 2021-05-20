# Trouble-shoot Erin Ice's Example for Smith et al's 2014 marginalized two-part
# model code.
#
# Date: July 18, 2019
# Author: James Henderson, PhD (jbhender@umich.edu)

# libraries: ------------------------------------------------------------------
library(tidyverse); library(foreign); 

# This library has a function for estimating log-normal means: ----------------
library(EnvStats)

# source code: ----------------------------------------------------------------
source('./marginal_mean_twopart.R')

# data for testing: -----------------------------------------------------------
atus_long <- read.dta("~/Dropbox/marginal_mean_twopart_ice/atus_longmerge_r.dta")
atus_long <- atus_long %>% filter(year >= 2005)

# create some normalized weights based on "wt06": -----------------------------
atus_long = mutate(atus_long,
                   weight = wt06 / {mean( wt06 ) + sum( wt06 - mean(wt06) ) } )
                   
# Fit the model as in the original function, using Nelder-Mead to optimize: ---
tm0 = system.time({
 fit1a = marginal_mean_twopart(
   totalcare ~ factor(year) + male*ba + age + raceeth + weekend | 
              factor(year) + male*ba + age + raceeth + weekend, 
              data = atus_long, method = 'BFGS',
              optim.control = list(maxit = 1e5))
 })
fit1a$counts
if( fit1a$convergence != 0 ) cat('Model did not converge.\n')

# Optimize using the BFGS algorithm and the (corrected) gradient: -------------
tm1 = system.time({
  fit1b = marginal_mean_twopart(
  totalcare ~ factor(year) + male*ba + age + raceeth + weekend | 
    factor(year) + male*ba + age + raceeth + weekend, 
  data = atus_long, method = 'BFGS',
#  weights = atus_long$wt06, ## Won't work as causes integer overflow
  optim.control = list(maxit = 1e5))
})

tm2 = system.time({
  fit1c = marginal_mean_twopart(
    totalcare ~ factor(year) + male*ba + age + raceeth + weekend | 
      factor(year) + male*ba + age + raceeth + weekend, 
    data = atus_long, method = 'BFGS',
    weights = atus_long$weight,
    optim.control = list(maxit = 1e5))
})

fit1b$counts
if( fit1b$convergence != 0 ) cat('Model did not converge.\n')

# You can compare tm0 and tm1 to see that the latter is much quicker.  

# Estimated Coefficients
bhat = fit1b$par
round( cbind(fit1a$par, fit1b$par), 3)
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

## Years only, let's check if this is returning sensible results: -------------
fit2a = marginal_mean_twopart(
  totalcare ~ 1 | 1, 
  data = atus_long, 
  weights = atus_long$weight,
  method = 'BFGS',
  optim.control = list(maxit = 1e5))

print(fit2a)
summary(fit2a)

mm_hat = predict( fit2a, atus_long, type = 'marginal_mean')
p_hat = predict( fit2a, atus_long, type = 'prob')
cm_hat = predict( fit2a, atus_long, type = 'conditional_mean')

#m2 = two_part(  totalcare ~ factor(year) | factor(year), 
#           data = atus_long)
#Y = m2$Y; X = m2$X; Z = m2$Z; w = 1

fit2b = marginal_mean_twopart(
  totalcare ~ 1 | 1 , 
  data = atus_long, 
  method = 'BFGS',
  optim.control = list(maxit = 1e5) )



## Summarize data at annual level
atus_annual = 
  atus_long %>% #group_by(year) %>% 
  summarize( p_nonzero = mean(totalcare > 0), 
             mean_nonzero = mean(totalcare[totalcare > 0]),
             est_nonzero = elnormAlt(totalcare[totalcare > 0], method = 'mle')$parameters[1]
             )

atus_annual_w = 
  atus_long %>% #group_by(year) %>% 
  summarize( p_nonzero = mean(weight*{totalcare > 0}), 
             mean_nonzero = mean({weight*totalcare}[totalcare > 0])
  )
plogis( fit2b$par[1]); plogis( fit2a$par[1])
#cbind( atus_annual$p_nonzero, plogis( fit2b$par[1] + c(0, fit2b$par[2:13]) ) )
#cbind( atus_annual_w$p_nonzero, plogis( fit2a$par[1] + c(0, fit2a$par[2:13]) ) )


## Get predictions for zeros
predict_nonzero = function(fit, data){
  Z = model.matrix(~factor(year), data = data)
  as.numeric( plogis( Z %*% fit$par[grep('^zero', names(fit$par))] ) )
}

atus_annual = atus_annual %>% 
  mutate( phat_nonzero2 = predict_nonzero( fit2, atus_annual),
          phat_nonzero2b = predict_nonzero( fit2b, atus_annual)
          )

## Get predictions for marginalized means, then nonzero means
predict_marginal_mean = function(fit, data){
  X = model.matrix(~factor(year), data = data)
  as.numeric( exp( X %*% fit$par[grep('^mean', names(fit$par))]) )
}

atus_annual = atus_annual %>% 
  mutate( est_marginal_mean = predict_marginal_mean(fit2, atus_annual),
          est_marginal_meanb = predict_marginal_mean(fit2b, atus_annual)
  ) %>%
  mutate( est_nonzero2 = est_marginal_mean / phat_nonzero2,
          est_nonzero2b = est_marginal_mean / phat_nonzero2b
          )

ggplot(atus_annual, aes(x = est_nonzero, y = est_nonzero2) ) +
  geom_point()

ggplot(atus_annual, aes(x = mean_nonzero, y = est_nonzero2) ) +
  geom_point()

ggplot(atus_annual, aes(x = mean_nonzero, y = est_nonzero) ) +
  geom_point()
