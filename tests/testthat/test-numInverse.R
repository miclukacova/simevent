library(survival)
library(testthat)

test_that("simEventData finds the inverse in the right way",{
  set.seed(857)
  # Generate data
  at_risk <- function(events) {
    return(c(
      1,1,                          # Always at risk for event 0 and 1
      as.numeric(events[3] < 1),    # Can experience event 2 once
      as.numeric(events[4] < 1)))   # Can experience event 3 once
  }
  beta <- matrix(rnorm(24,0,1), ncol = 4, nrow = 6)
  data_test <- simEventData(6000, beta = beta, eta = c(0.1,0.1,0.2,0.1))

  # Transform data into tstart tstop format
  data_int <- IntFormatData(data_test)

  # Fit models
  survfit_cens <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + N2 + N3, data = data_int)
  survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + N2 + N3, data = data_int)
  survfit_oper <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + A0 + N3, data = data_int[N2 == 0])
  survfit_cov <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0 + A0 + N2, data = data_int[N3 == 0])

  # Compare confidence intervals and true values
  expect_true(confint(survfit_cens, level = 0.99)[1,1] <= beta[1,1] & beta[1,1] <= confint(survfit_cens, level = 0.99)[1,2])
  expect_true(confint(survfit_cens, level = 0.99)[2,1] <= beta[2,1] & beta[2,1] <= confint(survfit_cens, level = 0.99)[2,2])
  expect_true(confint(survfit_cens, level = 0.99)[3,1] <= beta[5,1] & beta[5,1] <= confint(survfit_cens, level = 0.99)[3,2])
  expect_true(confint(survfit_cens, level = 0.99)[4,1] <= beta[6,1] & beta[6,1] <= confint(survfit_cens, level = 0.99)[4,2])
  expect_true(confint(survfit_death, level = 0.99)[1,1] <= beta[1,2] & beta[1,2] <= confint(survfit_death, level = 0.99)[1,2])
  expect_true(confint(survfit_death, level = 0.99)[2,1] <= beta[2,2] & beta[2,2] <= confint(survfit_death, level = 0.99)[2,2])
  expect_true(confint(survfit_death, level = 0.99)[3,1] <= beta[5,2] & beta[5,2] <= confint(survfit_death, level = 0.99)[3,2])
  expect_true(confint(survfit_death, level = 0.99)[4,1] <= beta[6,2] & beta[6,2] <= confint(survfit_death, level = 0.99)[4,2])
  expect_true(confint(survfit_oper, level = 0.99)[1,1] <= beta[1,3] & beta[1,3] <= confint(survfit_oper, level = 0.99)[1,2])
  expect_true(confint(survfit_oper, level = 0.99)[2,1] <= beta[2,3] & beta[2,3] <= confint(survfit_oper, level = 0.99)[2,2])
  expect_true(confint(survfit_oper, level = 0.99)[3,1] <= beta[6,3] & beta[6,3] <= confint(survfit_oper, level = 0.99)[3,2])
  expect_true(confint(survfit_cov, level = 0.99)[1,1] <= beta[1,4] & beta[1,4] <= confint(survfit_cov, level = 0.99)[1,2])
  expect_true(confint(survfit_cov, level = 0.99)[2,1] <= beta[2,4] & beta[2,4] <= confint(survfit_cov, level = 0.99)[2,2])
  expect_true(confint(survfit_cov, level = 0.99)[3,1] <= beta[5,4] & beta[5,4] <= confint(survfit_cov, level = 0.99)[3,2])
})
