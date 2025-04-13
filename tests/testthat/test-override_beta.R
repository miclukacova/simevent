library(survival)
library(testthat)

test_that("override_beta works", {
  set.seed(857)
  # Generate data
  beta <- matrix(rnorm(4*6), ncol = 4)
  at_risk <- function(events) {
    return(c(
      1,1,                          # Always at risk for event 0 and 1
      as.numeric(events[3] < 2),    # Can experience event 2 twice
      as.numeric(events[4] < 1)))   # Can experience event 3 once
  }
  data_test <- simEventData(1000, override_beta = list("N2 > 1" = c("N1" = 2)),
                            beta = beta, at_risk = at_risk)

  # Transform data into tstart tstop format
  data_int <- IntFormatData(data_test)

  # Fit models
  survfit1 <- coxph(Surv(tstart, tstop, Delta == 1) ~ L0 + A0 +
                      N2 + as.numeric(N2 > 1) + N3, data = data_int)

  # Compare confidence intervals and true values
  expect_true(confint(survfit1, level = 0.99)[1,1] <= beta[1,2] & beta[1,2] <= confint(survfit1, level = 0.99)[1,2])
  expect_true(confint(survfit1, level = 0.99)[2,1] <= beta[2,2] & beta[2,2] <= confint(survfit1, level = 0.99)[2,2])
  expect_true(confint(survfit1, level = 0.99)[3,1] <= beta[5,2] & beta[5,2] <= confint(survfit1, level = 0.99)[3,2])
  expect_true(confint(survfit1, level = 0.99)[4,1] <= 2 & 2 <= confint(survfit1, level = 0.99)[4,2])
  expect_true(confint(survfit1, level = 0.99)[5,1] <= beta[6,2] & beta[6,2] <= confint(survfit1, level = 0.99)[5,2])
})
