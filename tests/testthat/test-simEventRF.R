#library(survival)
#library(testthat)
#
#test_that("simEventCox simulates data in the right way",{
#
#  set.seed(373)
#  # The observed data
#  beta = matrix(c(1, 0,0.5,1), ncol = 2, nrow = 2)
#  data <- simSurvData(N = 10^4, beta = beta, cens = 0)
#
#  # Fit RF
#  RF_fit <- randomForestSRC::rfsrc(Surv(Time, Delta) ~ L0 + A0, data = data)
#
#  # Then simulate new data:
#  new_data <- simEventRF(10^5, RF_fit, L0_old = data$L0, A0_old = data$A0)
#
#  # Fit new Cox models
#  cox_death <- survival::coxph(survival::Surv(Time, N1 == 1) ~ L0 + A0, data = new_data)
#
#  # Compute confidence intervals
#  ci_death <- confint(cox_death)
#
#  # Compare confidence intervals and true values
#  expect_true(all(beta[,2] >= ci_death[, 1] & beta[,2] <= ci_death[, 2]))
#})
