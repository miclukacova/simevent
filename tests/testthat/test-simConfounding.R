library(survival)
library(testthat)

test_that("simConfounding simulates data in the right way",{
  set.seed(858)
  # Generate data
  data_test <- simConfounding(500)

  # Transform data into tstart tstop format
  data_int <- IntFormatData(data_test, N_cols = 5:6)

  # Fit models
  survfit_cens <- coxph(Surv(tstart, tstop, Delta == 0) ~ L0 + L + A, data = data_int)
  survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + L + A, data = data_int)
  survfit_oper <- coxph(Surv(tstart, tstop, Delta == 2) ~ L0  + L, data = data_int[A == 0])
  survfit_cov <- coxph(Surv(tstart, tstop, Delta == 3) ~ L0 + A, data = data_int[L == 0])

  # Compare confidence intervals and true values
  expect_true(confint(survfit_cens, level = 0.99)[1,1] <= 0 & 0 <= confint(survfit_cens, level = 0.99)[1,2])
  expect_true(confint(survfit_cens, level = 0.99)[2,1] <= 0 & 0 <= confint(survfit_cens, level = 0.99)[2,2])
  expect_true(confint(survfit_cens, level = 0.99)[3,1] <= 0 & 0 <= confint(survfit_cens, level = 0.99)[3,2])
  expect_true(confint(survfit_death, level = 0.99)[1,1] <= 1 & 1 <= confint(survfit_death, level = 0.99)[1,2])
  expect_true(confint(survfit_death, level = 0.99)[2,1] <= 1 & 1 <= confint(survfit_death, level = 0.99)[2,2])
  expect_true(confint(survfit_death, level = 0.99)[3,1] <= -1 & -1 <= confint(survfit_death, level = 0.99)[3,2])
  expect_true(confint(survfit_oper, level = 0.99)[1,1] <= 1 & 1 <= confint(survfit_oper, level = 0.99)[1,2])
  expect_true(confint(survfit_oper, level = 0.99)[2,1] <= 1 & 1 <= confint(survfit_oper, level = 0.99)[2,2])
  expect_true(confint(survfit_cov, level = 0.99)[1,1] <= 1 & 1 <= confint(survfit_cov, level = 0.99)[1,2])
  expect_true(confint(survfit_cov, level = 0.99)[2,1] <= -0.5 & -0.5 <= confint(survfit_cov, level = 0.99)[2,2])
})
