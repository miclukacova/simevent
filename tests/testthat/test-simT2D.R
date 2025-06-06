library(survival)
library(testthat)

test_that("simT2D simulates data in the right way",{
  set.seed(857)
  # Generate data
  data_test <- simT2D(1000)

  # Transform data into tstart tstop format
  data_int <- IntFormatData(data_test, N_cols = 6)

  # Fit models
  survfit_cens <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + L, data = data_int)
  survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + L, data = data_int)
  survfit_cov <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0 + A0, data = data_int[L == 0])

  # Compare confidence intervals and true values
  expect_true(confint(survfit_cens, level = 0.99)[1,1] <= 0 & 0 <= confint(survfit_cens, level = 0.99)[1,2])
  expect_true(confint(survfit_cens, level = 0.99)[2,1] <= 0 & 0 <= confint(survfit_cens, level = 0.99)[2,2])
  expect_true(confint(survfit_cens, level = 0.99)[3,1] <= 0 & 0 <= confint(survfit_cens, level = 0.99)[3,2])
  expect_true(confint(survfit_death, level = 0.99)[1,1] <= 1 & 1 <= confint(survfit_death, level = 0.99)[1,2])
  expect_true(confint(survfit_death, level = 0.99)[2,1] <= 0 & 0 <= confint(survfit_death, level = 0.99)[2,2])
  expect_true(confint(survfit_death, level = 0.99)[3,1] <= 1 & 1 <= confint(survfit_death, level = 0.99)[3,2])
  expect_true(confint(survfit_cov, level = 0.99)[1,1] <= 1 & 1 <= confint(survfit_cov, level = 0.99)[1,2])
  expect_true(confint(survfit_cov, level = 0.99)[2,1] <= 0 & 0 <= confint(survfit_cov, level = 0.99)[2,2])
})
