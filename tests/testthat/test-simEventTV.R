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
      val1 <- inverseScHazTV(p = p, t = t,
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
      val1 <- inverseScHazTV(p = p, t = t,
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

test_that("That the simEventTV simulates correctly when no time varrying effects",{
  set.seed(34730)

  eta <- rep(0.1, 2)
  term_deltas <- c(0)
  tv_eff <- matrix(0.1, ncol = 2, nrow = 4)
  beta <- matrix(nrow = 4, ncol = 2, 0.1)
  t_prime <- 1000

  data <- simEventTV(N = 5*10^3, t_prime = t_prime, tv_eff = tv_eff, eta = eta,
                     term_deltas = term_deltas, beta = beta, lower = 10^(-20), upper = 10^2,
                     max_events = 5)

  data <- IntFormatData(data, N_cols = 6:7)

  survfit0 <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + N1, data = data)
  survfit1 <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + N1, data = data)

  errors <- sum(abs(survfit0$coefficients - beta[c(1,2,4),1]) +
    abs(survfit1$coefficients - beta[c(1,2,4),2]))

  t_prime <- 1
  tv_eff <- matrix(0.1, ncol = 2, nrow = 4)

  data <- simEventTV(N = 2*10^4, t_prime = t_prime, tv_eff = tv_eff, eta = eta,
                     term_deltas = term_deltas, beta = beta, lower = 10^(-20), upper = 10^2,
                     max_events = 5)

  data <- IntFormatData(data, N_cols = 6:7)

  survfit0 <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + N1, data = data)
  survfit1 <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + N1, data = data)

  errors <- errors + sum(abs(survfit0$coefficients - beta[c(1,2,4),1]) +
    abs(survfit1$coefficients - beta[c(1,2,4),2]))

  expect_equal(errors, 0, tolerance = 0.1*12)
})

test_that("That the simEventTV simulates correctly when time varrying effects",{
  set.seed(623547)
  eta <- rep(0.1, 2)
  term_deltas <- c(0)
  beta <- matrix(nrow = 4, ncol = 2, 0.1)
  tv_eff <- matrix(0.3, ncol = 2, nrow = 4)
  t_prime <- 1

  data <- simEventTV(N = 5*10^3, t_prime = t_prime, tv_eff = tv_eff, eta = eta,
                     term_deltas = term_deltas, beta = beta, lower = 10^(-20), upper = 10^2,
                     max_events = 5)

  data <- IntFormatData(data, timeVar = TRUE, t_prime = t_prime, N_cols = 6:7)

  survfit0 <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0:strata(t_group) + A0:strata(t_group)
                    + N1:strata(t_group), data = data)
  survfit1 <- coxph(Surv(tstart, tstop, Delta == 1)~ L0:strata(t_group) + A0:strata(t_group)
                    + N1:strata(t_group), data = data)

  errors <- sum(abs(survfit0$coefficients[c(1,3,5)] - beta[c(1,2,4),1]))+
    sum(abs(survfit1$coefficients[c(1,3,5)] - beta[c(1,2,4),2]))

  betap <- beta + tv_eff
  errors <- errors + sum(abs(survfit0$coefficients[c(2,4,6)] - betap[c(1,2,4),1]))+
    sum(abs(survfit1$coefficients[c(2,4,6)] - betap[c(1,2,4),2]))
  expect_equal(errors, 0, tolerance = 0.1*12)
})

