#include "Rcpp.h"
#include <Rinternals.h>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector cpp_ripcliques_to_jtree_zo(Function print, NumericMatrix cliques, NumericMatrix jtree, int p, int num_cliques )
{
  NumericVector score(p);
  if ( num_cliques >= 2)
  {
    for ( int i = 2 ; i < num_cliques + 1 ; i ++)
    {
      int sum;
      int max = 0;
      int index = 0;
      for (int k = 0 ; k < (i-1) ; k++)
      {
        sum = 0;
        for (int aux = 0 ; aux < p ; aux ++)
        {
          if(cliques.column(i-1)[aux]!=0)
          {
            if(cliques.column(k)[aux]!=0)
            {
              sum = sum + cliques.column(i-1)[aux]*cliques.column(k)[aux];
            }
          }
        }
        score[k] = sum;
        if (sum > max)
        {
          max = sum;
          index = k;
        }
      }
      if (max != 0)
      {
        jtree.row(i-1)[index] = 1;
        jtree.row(index)[i-1] = 1;
        
      }
    }
  }
  
  return jtree;
}

