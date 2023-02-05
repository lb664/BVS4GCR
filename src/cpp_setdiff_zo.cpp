#include "Rcpp.h"
#include <Rinternals.h>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector cpp_setdiff_zo(NumericVector a , NumericVector b )
{
  if(b.size() > 0)
  {
    for (int j = 0 ; j < b.size() ; j++)
    {
      if (b[j] < 0)
      {
        b[j] = a[j];
      }
    }
  }
  return (b);
}
