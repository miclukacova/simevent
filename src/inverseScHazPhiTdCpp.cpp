#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double inverseScHazPhiTdCpp(
    double p,
    double t,
    NumericVector T_star,
    double lower,
    double upper,
    NumericVector eta,
    NumericVector nu,
    NumericVector beta2,
    NumericVector phi0,
    NumericVector at_risk,
    double tol = 1e-9,
    int max_iter = 100) {

  // Closed-form cumulative hazard
  auto cum_haz = [&](double u) {

    double total = 0.0;

    for (int k = 0; k < eta.size(); ++k) {

      if (at_risk[k] == 0)
        continue;

      double shape = nu[k];
      double rate  = beta2[k];

      double a = rate * t;
      double b = rate * (t + u);

      // Lower incomplete gamma via pgamma
      double G_b =
        R::pgamma(b, shape, 1.0, 1, 0);

      double G_a =
        R::pgamma(a, shape, 1.0, 1, 0);

      double integral;

      if(std::abs(rate) < 1e-12) {

        // Standard Weibull cumulative hazard if beta_2 = 0
        integral =
          (std::pow(t + u, shape) -
          std::pow(t, shape)) / shape;

      } else {

        double a = rate * t;
        double b = rate * (t + u);

        double G_b =
          R::pgamma(b, shape, 1.0, 1, 0);

        double G_a =
          R::pgamma(a, shape, 1.0, 1, 0);

        integral =
          std::pow(rate, -shape) *
          R::gammafn(shape) *
          (G_b - G_a);
      }

      double const_part =
        eta[k] *
        shape *
        phi0[k] *
        std::exp(rate * T_star[k]);

      total +=
        at_risk[k] *
        const_part *
        integral;
    }

    return total;
  };

  // Root finding via bisection
  double a = lower;
  double b = upper;

  double fa = cum_haz(a) - p;
  double fb = cum_haz(b) - p;

  if (fa * fb > 0) {

    Rcout << "cum_haz(a) = " << cum_haz(a)
          << ", cum_haz(b) = " << cum_haz(b)
          << ", p = " << p
          << std::endl;

    stop("Root not bracketed");
  }

  for (int iter = 0; iter < max_iter; ++iter) {

    double mid = 0.5 * (a + b);
    double fmid = cum_haz(mid) - p;

    if (std::abs(fmid) < tol)
      return mid;

    if (fa * fmid < 0) {
      b = mid;
      fb = fmid;
    } else {
      a = mid;
      fa = fmid;
    }
  }

  warning("Maximum iterations reached");

  return 0.5 * (a + b);
}
