#include "Rcpp.h"
#include <Rinternals.h>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::export]]

List cpp_check_chordal(NumericMatrix g, NumericVector capU, int p, int chordal, NumericVector order)
{
  for (int i = 1 ; i < p ; i ++)
  {
    NumericVector score(p-i);
    int max = 0;
    int index_max = 0;
    
    for (int ui = 0 ; ui < (p-i) ; ui++)
    {
      int iterator  = 0;
      int u = capU[ui]-1;
      for(int k = 0 ; k < g.row(u).size() ; k++)
      {
        if (g.row(u)[k]!=0)
        {
          bool boolean = TRUE;
          int l = 0;
          while (l < (p-i) && boolean)
          {
              if (k+1==capU(l))
              {
                  boolean = FALSE;
              }
              l ++;
          }
          if (boolean)
          {
              iterator  ++;
          }
        }
      }
      if(iterator > max)
      {
        max = iterator ;
        index_max = ui;
      }
      score[ui] = iterator ;
    }
    order[i]=capU(index_max);
    capU.erase(index_max);
    bool chordal_1_or_0 = TRUE;
    NumericVector pa(i+1);
    int index_pa = 0;
    for(int k = 0 ; k < g.row(order[i]-1).size() ; k++)
    {
      if (g.row(order[i]-1)[k]!=0)
      {
        int z = 0 ;
        while (z<i)
        {
          if (order[z]==k+1)
          {
            pa[index_pa] = k;
            index_pa++;
            if (g.row(k)[k]!=1)
            {
              chordal_1_or_0 = FALSE;
              break;
            }
            if (index_pa-1>0)
            {
              for (int l = 0 ; l < index_pa-1 ; l++)
              {
                if (g.row(k)[pa[l]]!=1 || g.row(pa[l])[k]!=1)
                {
                  chordal_1_or_0 = FALSE;
                  break;
                }
              }
            }
            if (!chordal_1_or_0)
            {
              break;
            }
          }
          
          z++;
        }
        if (!chordal_1_or_0)
        {
          break;
        }
      }
    }
    
    if (!chordal_1_or_0)
    {
      chordal = 0;
      break;
    }
  }
  List output = List::create(chordal,order);
  return (output) ;
}
