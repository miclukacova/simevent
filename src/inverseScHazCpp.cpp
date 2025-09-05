#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double inverseScHazCpp(double p,
                        double t,
                        double lower,
                        double upper,
                        NumericVector eta,
                        NumericVector nu,
                        NumericVector phi,
                        NumericVector at_risk,
                        double tol = 1e-9,
                        int max_iter = 100) {

  auto cum_haz = [&](double u) {
    double sum = 0.0;
    for (int k = 0; k < eta.size(); ++k) {
      sum += at_risk[k] * eta[k] * phi[k] * (pow(t + u, nu[k]) - pow(t, nu[k]));
    }
    return sum;
  };

  double a = lower;
  double b = upper;
  double fa = cum_haz(a) - p;
  double fb = cum_haz(b) - p;

  if (fa * fb > 0) {
    stop("Function does not bracket root: adjust upper and lower");
  }

  for (int iter = 0; iter < max_iter; ++iter) {
    double mid = 0.5 * (a + b);
    double fmid = cum_haz(mid) - p;

    if (std::abs(fmid) < tol) {
      return mid;
    }

    if (fa * fmid < 0) {
      b = mid;
      fb = fmid;
    } else {
      a = mid;
      fa = fmid;
    }
  }

  warning("Max iterations reached in inverseCumHazard");
  return 0.5 * (a + b);
}
