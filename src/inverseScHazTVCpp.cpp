#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double inverseScHazTVCpp(double p,
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

  // Cumulative hazard function with time-varying components
  auto cum_haz = [&](double u) {
    double sum = 0.0;
    if(u < t_prime) {
      // Before t_prime the cumulative hazard is computed in the usual way
      for (int k = 0; k < eta.size(); ++k) {
        sum += at_risk[k] * eta[k] * phi[k] * pow(u, nu[k]);
      }
    }
    else{
      // After t_prime the cumulative hazard changes due to a change in phi
      for (int k = 0; k < eta.size(); ++k) {
        sum += at_risk[k] * eta[k] * (phi[k] - phi_prime[k]) * pow(t_prime, nu[k]) +
          at_risk[k] * eta[k] * phi_prime[k] * pow(u, nu[k]);
      }
    }
    return sum;
  };

  // Cumulative hazard of weighting time
  auto cum_haz_wait = [&](double u, double t) {
    return cum_haz(u + t) - cum_haz(t);
  };

  double a = lower;
  double b = upper;
  double fa = cum_haz_wait(a, t) - p;
  double fb = cum_haz_wait(b, t) - p;

  // Check root is bracketed
  if (fa * fb > 0) {
    stop("Function does not bracket root: adjust upper and lower");
  }

  // Bisection method to find root of the cumulative hazard of the weighting time
  for (int iter = 0; iter < max_iter; ++iter) {
    double mid = 0.5 * (a + b);
    double fmid = cum_haz_wait(mid, t) - p;

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

  warning("Max iterations reached in inverseScHazTVCpp");
  return 0.5 * (a + b);
}
