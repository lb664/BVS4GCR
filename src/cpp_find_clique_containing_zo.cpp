#include "Rcpp.h"
#include <Rinternals.h>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::export]]

int cpp_find_clique_containing_zo(int a, NumericMatrix cliques)
{
  int output = -1;
  NumericVector aux = cliques.row(a-1);
  for (int i = 0 ; i < aux.size() ; i++)
  {
    if (aux[i]!=0)
    {
      output = i + 1;
      break;
    }
  }
  return output;
}
