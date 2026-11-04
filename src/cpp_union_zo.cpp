#include "Rcpp.h"
#include <Rinternals.h>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector cpp_union_zo(NumericVector b)
{
  if(b.size() > 0)
  {
    for (int j = 0 ; j < b.size() ; j++)
    {
      if (b[j] > 1)
      {
        b[j] = 1;
      }
    }
  }
  return (b);
}
