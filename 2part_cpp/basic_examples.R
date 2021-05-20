# Hello World Examples

# libraries: ------------------------------------------------------------------
library(Rcpp)

# The basic example: ----------------------------------------------------------
sourceCpp('./timesTwo.cpp')
timesTwo(1:5)

# Adding console messages with Rcout (useful for debugging): ------------------

sourceCpp(code =
'
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector timesThree(NumericVector x) {
  Rcout << "Multiplying by 3:"; 
  return x * 3;
}
'
)
timesThree(1:5)

# illustrating use of namespace: ----------------------------------------------
sourceCpp(code =
'
#include <Rcpp.h>
// [[Rcpp::export]]
Rcpp::NumericVector timesFour(Rcpp::NumericVector x) {
  Rcpp::Rcout << "Multiplying by 4:"; 
  return x * 4;
}
'
)
timesFour(1:5)
