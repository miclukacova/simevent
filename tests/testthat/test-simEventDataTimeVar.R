library(testthat)
library(survival)

test_that("That inverseScHaz and inverseScHazTimeVar are qual when no time varying effects",{
  eta <- c(0.1,0.1)
  nu <- c(1.1,1.1)
  phi <- c(1,1)
  at_risk <- c(1, 1)
  phi_prime <- c(2, 2)
  t_prime <- 10000

  t_vals <- seq(0.1, 10, length.out = 5)
  p_vals <- seq(0.1, 0.9, length.out = 5)

  for (t in t_vals) {
    for (p in p_vals) {
      val1 <- inverseScHazTimeVar(p = p, t = t,
                                  eta = eta, nu = nu, phi = phi,
                                  at_risk = at_risk,
                                  lower = 1e-15, upper = 200,
                                  t_prime = t_prime,
                                  phi_prime = phi_prime)

      val2 <- inverseScHaz(p = p, t = t,
                           eta = eta, nu = nu, phi = phi,
                           at_risk = at_risk,
                           lower = 1e-15, upper = 200)

      msg <- paste("Mismatch at t =", t, "and p =", p)
      expect_equal(val1, val2, tolerance = 1e-6, info = msg)
    }
  }

  t_prime <- 0.1
  phi_prime <- c(1,1)

  for (t in t_vals) {
    for (p in p_vals) {
      val1 <- inverseScHazTimeVar(p = p, t = t,
                                  eta = eta, nu = nu, phi = phi,
                                  at_risk = at_risk,
                                  lower = 1e-15, upper = 200,
                                  t_prime = t_prime,
                                  phi_prime = phi_prime)

      val2 <- inverseScHaz(p = p, t = t,
                           eta = eta, nu = nu, phi = phi,
                           at_risk = at_risk,
                           lower = 1e-15, upper = 200)

      msg <- paste("Mismatch at t =", t, "and p =", p)
      expect_equal(val1, val2, tolerance = 1e-6, info = msg)
    }
  }
})

#test_that("That the simEventDataTimeVar simulates correctly",{
#  set.seed(34729)
#  eta <- rep(0.1, 4)
#  term_deltas <- c(0,1)
#  time_var_eff <- matrix(0.1, ncol = 4, nrow = 6)
#  beta <- matrix(nrow = 6, ncol = 4, rnorm(4*6, sd = 0.5))
#  t_prime <- 1
#  time_var_eff <- matrix(0.3, ncol = 4, nrow = 6)
#
#  data <- simEventDataTimeVar(N = 3*10^3, t_prime = t_prime, time_var_eff = time_var_eff, eta = eta,
#                              term_deltas = term_deltas, beta = beta, lower = 10^(-20), upper = 10^2,
#                              max_events = 5)
#
#  data <- IntFormatData(data, time_var = TRUE, t_prime = t_prime)
#
#  survfit0 <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0:strata(t_group) + A0:strata(t_group)
#                    + N2:strata(t_group) + N3:strata(t_group), data = data_new)
#  survfit1 <- coxph(Surv(tstart, tstop, Delta == 1)~ L0:strata(t_group) + A0:strata(t_group)
#                    + N2:strata(t_group) + N3:strata(t_group), data = data_new)
#  survfit2 <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0:strata(t_group) + A0:strata(t_group)
#                    + N2:strata(t_group) + N3:strata(t_group), data = data_new)
#  survfit3 <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0:strata(t_group) + A0:strata(t_group)
#                    + N2:strata(t_group) + N3:strata(t_group), data = data_new)
#
#
#  any(survfit0$coefficients[c(1,3,5,7)] - beta[c(1,2,5,6),1] > 0.2)
#  any(survfit1$coefficients[c(1,3,5,7)] - beta[c(1,2,5,6),2] > 0.2)
#  any(survfit2$coefficients[c(1,3,5,7)] - beta[c(1,2,5,6),3] > 0.2)
#  any(survfit3$coefficients[c(1,3,5,7)] - beta[c(1,2,5,6),4] > 0.2)
#
#  betap <- beta + time_var_eff
#  any(survfit0$coefficients[c(2,4,6,8)] - betap[c(1,2,5,6),1] > 0.2)
#  any(survfit1$coefficients[c(2,4,6,8)] - betap[c(1,2,5,6),2] > 0.2)
#  any(survfit2$coefficients[c(2,4,6,8)] - betap[c(1,2,5,6),3] > 0.2)
#  any(survfit3$coefficients[c(2,4,6,8)] - betap[c(1,2,5,6),4] > 0.2)
#
#})

