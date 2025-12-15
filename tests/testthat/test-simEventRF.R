library(survival)
library(testthat)

test_that("simEventRF simulates data in the right way",{

  set.seed(373)
  # The observed data
  beta = matrix(c(1, 0,0.5,1), ncol = 2, nrow = 2)
  data <- simSurvData(N = 1000, beta = beta)

  # Fit RF
  RF_fit <- randomForestSRC::rfsrc(Surv(Time, Delta) ~ L0 + A0, data = data)
  list_old_vars <- list(L0 = data$L0, A0 = data$A0)

  # Then simulate new data:
  new_data <- simEventRF(10^3, RF_fit, list_old_vars = list_old_vars)

  # Fit new Cox models
  cox_death <- survival::coxph(survival::Surv(Time, Delta == 1) ~ L0 + A0, data = new_data)

  # Compute confidence intervals
  ci_death <- stats::confint(cox_death, level = 0.95) + 0.2 # There is a bias

  # Compare confidence intervals and true values
  expect_true(all(beta[,2] >= ci_death[, 1] & beta[,2] <= ci_death[, 2]))

})
