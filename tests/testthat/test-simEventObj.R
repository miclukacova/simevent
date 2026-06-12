library(survival)
library(testthat)

test_that("simEventObj simulates data in the right way",{

  print("We are skipping this check, fix this at some point")

  if(FALSE){
    set.seed(926489)
    # The observed data
    beta = matrix(c(0.5,-1,-0.5,0.5,0,0.5), ncol = 3, nrow = 2)
    data <- simCRdata(N = 10^4, beta = beta)
    old_vars <- data[, c("L0", "A0")]

    # Fit Cox Models
    cox1 <- survival::coxph(survival::Surv(Time, Delta == 1) ~ L0 + A0, data = data)
    cox2 <- survival::coxph(survival::Surv(Time, Delta == 2)~ L0 + A0, data = data)

    # Create Object
    cox_fits <- list("D" = cox1, "L" = cox2)
    class(cox_fits) <- "simevent"

    # Equip with predict2 method
    predict2 <- function(obj, ...) {
      UseMethod("predict2")
    }
    predict2.simevent <- function(obj, sim_data){
      # Base hazards
      basehazz_list <- lapply(obj, function(model) basehaz(model, centered = FALSE))

      # Individual specific term
      cox_term <- lapply(obj, function(model)
        exp(stats::predict(model, newdata = sim_data, type="lp", reference = "zero")))

      # Chf
      chf_list <- lapply(1:length(obj), function(j) cox_term[[j]] %*% t(basehazz_list[[j]][["hazard"]]))
      chf <- array(dim = c(c(dim(chf_list[[1]])), length(obj)))
      for(j in 1:length(obj)) chf[,,j] <- chf_list[[j]]

      # Return list
      list(time = basehazz_list[[1]][["time"]], chf = chf)
    }

    # Simulate new data
    new_data <- simEventObj(10^4, cox_fits, old_vars = old_vars)

    # Check whether new data corresponds to the old data
    coxfit1 <- coxph(Surv(Time, Delta == 1) ~ L0 + A0, data = new_data)
    coxfit2 <- coxph(Surv(Time, Delta == 2) ~ L0 + A0, data = new_data)

    # Compute confidence intervals
    ci1 <- confint(cox1)
    coef1 <- coxfit1$coefficients

    ci2 <- confint(cox2)
    coef2 <- coxfit2$coefficients

    # Compare confidence intervals and true values
    expect_true(all(coef1 >= ci1[, 1] & coef1 <= ci1[, 2]))
    expect_true(all(coef2 >= ci2[, 1] & coef2 <= ci2[, 2]))
  }

})
