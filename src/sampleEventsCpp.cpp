#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector sampleEventsCpp(NumericMatrix probs) {
  int n = probs.ncol();
  int k = probs.nrow();
  IntegerVector result(n);

  for (int j = 0; j < n; ++j) {
    NumericVector p(k);
    double total = 0.0;

    for (int i = 0; i < k; ++i) {
      p[i] = probs(i, j);
      total += p[i];
    }

    for (int i = 0; i < k; ++i) {
      p[i] /= total;
    }

    double u = R::runif(0.0, 1.0);
    double cumulative = 0.0;

    for (int i = 0; i < k; ++i) {
      cumulative += p[i];
      if (u <= cumulative) {
        result[j] = i;
        break;
      }
    }
  }

  return result;
}
