library(survival)
library(testthat)

test_that("simEventCox simulates data in the right way",{

  set.seed(92648)
  # The observed data
  data_obs <- simDisease(N = 5000)
  data_obs <- IntFormatData(data_obs, N_cols = 6)

  # Fit some Cox models
  cox_death <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + L, data = data_obs)
  cox_Disease <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + A0, data = data_obs[L == 0])

  # Then simulate new data:
  cox_fits <- list("D" = cox_death, "L" = cox_Disease)
  new_data <- simEventCox(1000, cox_fits = cox_fits, L0_old = data_obs$L0, A0_old = data_obs$A0)
  new_data <- IntFormatData(new_data, N_cols = 6:7)

  # Fit new Cox models
  cox_death2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + L, data = new_data)
  cox_Disease2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + A0, data = new_data[L == 0])

  # Compute confidence intervals
  ci_death <- confint(cox_death2)
  coef_death <- cox_death$coefficients

  ci_Disease <- confint(cox_Disease2)
  coef_Disease <- cox_Disease$coefficients

  # Compare confidence intervals and true values
  expect_true(all(coef_death >= ci_death[, 1] & coef_death <= ci_death[, 2]))
  expect_true(all(coef_death >= ci_death[, 1] & coef_death <= ci_death[, 2]))
})
