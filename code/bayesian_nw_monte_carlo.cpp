#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// Monte Carlo simulation of species extinctions
// Input:
// - A: Adjacency matrix of the food web
// - Pb: Vector of baseline extinction probabilities
// - alpha: First parameter of the Beta distribution
// - beta: Second parameter of the Beta distribution
// - nreps: Number of iterations for Monte Carlo sim of Bayesian network
// Output:
// - A logical matrix, with columns corresponding to the species and rows to
//   replicate Monte Carlo simulations. The (r,i)th entry of the matrix is TRUE
//   if species i survived in replicate r; otherwise it is FALSE.
// [[Rcpp::export]]
LogicalMatrix persistence_mc(NumericMatrix A, NumericVector Pb,
                             double alpha, double beta, int nreps) {
  int S = Pb.size(), ind_basal = 0, i, j, rep;
  double A_times_extant, frac, p_ext;
  LogicalMatrix extant(nreps,S);
  NumericVector rowsums(S), rand(S);
  for (i = 0; i < S; i++) {
    rowsums(i) = sum(A(i,_));
    if (rowsums(i) == 0.0) ind_basal++;
  }
  for (rep = 0; rep < nreps; rep++) {
    rand = Rcpp::runif(S);
    for (i = 0; i < ind_basal; i++) extant(rep,i) = (rand(i) < Pb(i) ? 0.0 : 1.0);
    for (i = ind_basal; i < S; i++) {
      A_times_extant = 0.0;
      for (j = 0; j < S; j++) A_times_extant += A(i,j) * extant(rep,j);
      frac = 1.0 - A_times_extant / rowsums(i);
      p_ext = Pb(i) + (1 - Pb(i)) * R::pbeta(frac, alpha, beta, true, false);
      extant(rep,i) = (rand(i) < p_ext ? 0.0 : 1.0);
    }
  }
  return(extant);
}
