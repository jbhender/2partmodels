#include <Rcpp.h>
#include <math.h>  // exp, log
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double mm2p_loglik_cpp(NumericVector & y,  NumericVector & w,
                       NumericVector & p0, NumericVector & v,
                       double & s) {
  double out = 0;
  int npos = 0; 
  static const double pi = 3.141593; 
  const double v1 = .5 * log(2 * pi) + s; 
  const double v2 = exp( 2 * s );
  double v3; 
  
  for (int i = 0; i < p0.size(); i++) 
    {
      out += -w[i] * log( 1 + exp(p0[i]));
      if ( y[i] > 0 )
      {
        npos++;
        v3 = log(y[i]) + p0[i] - log( 1 + exp(p0[i]) ) + .5 * v2 - v[i]; 
        out += w[i] * ( p0[i] - log(y[i]) - v1 - .5 / v2 * v3 * v3 ); 
      }
    }

  return out;
}

/* Orignal R code
 loglik = 
   sum( -w * log( 1 + exp(p0) ) ) +
   sum( w[pos] * { p0[pos] - log( y[pos] ) - 
                   .5 * log( 2*pi ) - s - #v1
                   .5 / exp(2 * s) * #v2
                   { log( y[pos] ) + p0[pos] - log( 1 + exp(p0[pos]) ) +
                     .5 * exp(2 * s) - #v2
                     v[pos]  
                   }^2 
                 } 
  )
*/

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
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

# key values for the likelihood computation below
p0 = Z %*% coef(zero)
v = X %*% coef(cont)
s = sigma(cont)
y = Y
    
# log likelihood
pos = which(y > 0)
loglik = sum( -w * log( 1 + exp(p0) ) ) +
  sum( w[pos]*{p0[pos] - log( y[pos] ) - .5*log( 2*pi ) - s -
      .5 / {exp(s)^2 } * { log( y[pos] ) + p0[pos] -
          log( 1 + exp(p0[pos]) ) +
          exp(s)^2*.5 - v[pos]  }^2} 
  )
loglik_cpp = mm2p_loglik_cpp(y, w, p0, v, s)
c( 'loglik' = loglik, 'loglik_cpp' = loglik_cpp)
stopifnot( abs(loglik - loglik) < 1e-4 )

## Timing comparisons
y = rep(y, 100); w = rep(w, 100); p0 = rep(p0, 100); v = rep(v, 100)
loglik_R = function(y, w, p0, v, s) {
  pos = which(y > 0)
  sum( -w * log( 1 + exp(p0) ) ) +
    sum( w[pos]*{p0[pos] - log( y[pos] ) - .5*log( 2*pi ) - s -
                 .5 / {exp(s)^2 } * 
                 { log( y[pos] ) + p0[pos] - log( 1 + exp(p0[pos]) ) + 
                   exp(s)^2*.5 - v[pos]  }^2 
                  } 
    )
  loglik
}

microbenchmark::microbenchmark(loglik_R)
microbenchmark::microbenchmark(mm2p_loglik_cpp(y, w, p0, v, s))
*/
