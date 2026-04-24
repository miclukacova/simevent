library(survival)
library(testthat)

test_that("simEventCox simulates data in the right way",{

  set.seed(92648)
  # The observed data
  data_obs <- simDisease(N = 5000)
  data_obs <- IntFormatData(data_obs, N_cols = 6)

  # Fit some Cox models
  cox_cens <- survival::coxph(survival::Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + L, data = data_obs)
  cox_death <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + L, data = data_obs)
  cox_Disease <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + A0, data = data_obs[L == 0])

  # Then simulate new data:
  cox_fits <- list("C" = cox_cens, "D" = cox_death, "L" = cox_Disease)
  list_old_vars <- list("L0" = data_obs$L0, "A0" = data_obs$A0)
  new_data <- simEventCox(500, cox_fits = cox_fits, list_old_vars = list_old_vars,
                          term_events = c(1,2), n_event_max = c(1,1,1))
  new_data <- IntFormatData(new_data, N_cols = 6:8)

  # Fit new Cox models
  cox_cens2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 0) ~ L0 + A0 + L, data = new_data)
  cox_death2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 1) ~ L0 + A0 + L, data = new_data)
  cox_Disease2 <- survival::coxph(survival::Surv(tstart, tstop, Delta == 2) ~ L0 + A0, data = new_data[L == 0])

  # Compute confidence intervals
  ci_cens <- confint(cox_cens2)
  coef_cens <- cox_cens$coefficients

  ci_death <- confint(cox_death2)
  coef_death <- cox_death$coefficients

  ci_Disease <- confint(cox_Disease2)
  coef_Disease <- cox_Disease$coefficients

  # Compare confidence intervals and true values
  expect_true(all(coef_cens >= ci_cens[, 1] & coef_cens <= ci_cens[, 2]))
  expect_true(all(coef_death >= ci_death[, 1] & coef_death <= ci_death[, 2]))
  expect_true(all(coef_Disease >= ci_Disease[, 1] & coef_Disease <= ci_Disease[, 2]))
})
