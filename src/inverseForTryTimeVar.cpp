#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double inverseScHazCppTryTimeVar(double p,
                                 double t,
                                 double lower,
                                 double upper,
                                 double t_prime,
                                 NumericVector eta,
                                 NumericVector nu,
                                 NumericVector phi,
                                 NumericVector phi_prime,
                                 NumericVector at_risk,
                                 double tol = 1e-9,
                                 int max_iter = 100) {

  auto cum_haz = [&](double u) {
    double sum = 0.0;
    if(u < t_prime) {
      for (int k = 0; k < eta.size(); ++k) {
        sum += at_risk[k] * eta[k] * phi[k] * pow(u, nu[k]);
      }
    }
    else{
      for (int k = 0; k < eta.size(); ++k) {
        sum += at_risk[k] * eta[k] * (phi[k] - phi_prime[k]) * pow(t_prime, nu[k]) +
          at_risk[k] * eta[k] * phi_prime[k] * pow(u, nu[k]);
      }
    }
    return sum;
  };

  auto cum_haz_wait = [&](double u, double t) {
    return cum_haz(u + t) - cum_haz(t);
  };

  double a = lower;
  double b = upper;
  double fa = cum_haz_wait(a) - p;
  double fb = cum_haz_wait(b) - p;

  if (fa * fb > 0) {
    stop("Function does not bracket root: adjust upper and lower");
  }

  for (int iter = 0; iter < max_iter; ++iter) {
    double mid = 0.5 * (a + b);
    double fmid = cum_haz_wait(mid) - p;

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
