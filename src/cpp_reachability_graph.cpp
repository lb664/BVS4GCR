#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::vec find_and_replace(arma::vec A,
                           double find_val = 3.0,
                           double replace_val = 1.5) {
  arma::uvec idx = find(A > find_val); // Substitute == with >, <, >=, <=, !=
  A.elem(idx).fill(replace_val);        // Retrieve elements from positional index
  // Fill with the appropriate value
  return A;
}

// [[Rcpp::export]]
arma::mat cpp_reachability_graph(arma::mat g)
{
  int p = g.n_rows;
  arma::mat A = g;
  arma::mat reach_graph = arma::zeros(p,p);
  int i;
  
  for(i=0; i<(p-1); i++){
    reach_graph = reach_graph + A;
    A = A*g;
  }
  
  return reach_graph;
}


