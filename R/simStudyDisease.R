#' Simulation Study of the Kaplan-Meier Estimator Post Disease Diagnosis
#'
#' This function performs a simulation study evaluating the Kaplan-Meier estimator
#' in a disease setting Data is generated using `simDisease`, with
#' right censoring incorporated. The covariate `L1` is omitted. The estimand is:
#'
#' \deqn{S(\tau | A_0 = a_0) = P(T \geq \tau | A_0 = a_0), \quad a_0 \in \{0, 1\}, \ \tau \in \mathbb{R}}
#'
#' Here, \eqn{T} is the survival time post Disease diagnosis. After Disease onset, only death and censoring occur (no competing risks),
#' making the Kaplan-Meier estimator unbiased for the survival probability.
#'
#' For each of the \code{B} simulation runs, a dataset with \code{N} individuals is generated. The Kaplan-Meier estimate,
#' standard error (SE), and 95% confidence interval (CI) for the survival function at time \code{tau} are computed
#' separately for \eqn{A_0 = 0} and \eqn{A_0 = 1}.
#'
#' @title simStudyDisease
#'
#' @param N Integer. Number of individuals in each simulated dataset.
#' @param B Integer. Number of simulation repetitions.
#' @param beta_L0_D Numeric. Effect of L0 on risk of death.
#' @param beta_L0_L Numeric. Effect of L0 on risk of Disease.
#' @param beta_A0_L Numeric. Effect of A0 on risk of Disease.
#' @param beta_L_D Numeric. Effect of time-varying covariate on death risk.
#' @param beta_A0_D Numeric. Effect of A0 on risk of death.
#' @param eta Numeric vector of length 3. Shape parameters for the Weibull hazard \eqn{\eta \nu t^{\nu - 1}}.
#'   Default is \code{rep(0.1, 3)}.
#' @param nu Numeric vector of length 3. Scale parameters for the Weibull hazard. Default is \code{rep(1.1, 3)}.
#' @param tau Numeric. Time point at which the survival probability is estimated.
#'
#' @return A matrix of dimension \code{B x 8}, with the following columns:
#' \item{Est 0:}{ KM estimate of survival at \code{tau} for A_0 = 0}
#' \item{Est 1:}{ KM estimate of survival at \code{tau} for A_0 = 1}
#' \item{SE 0:}{ Standard error of the KM estimate for A_0 = 0}
#' \item{SE 1:}{ Standard error of the KM estimate for A_0 = 1}
#' \item{Lower 0:}{ Lower 95 % CI bound for A_0 = 0}
#' \item{Lower 1:}{ Lower 95 % CI bound for A_0 = 1}
#' \item{Upper 0:}{ Upper 95 % CI bound for A_0 = 0}
#' \item{Upper 1:}{ Upper 95 % CI bound for A_0 = 1}
#'
#' @export
#'
#' @examples
#' simStudyDisease(N = 100, B = 10, beta_L0_D = 1, beta_L0_L = 1, beta_A0_L = -1,
#'             beta_L_D = 1, beta_A0_D = 0, tau = 1)
simStudyDisease <- function(N, B, beta_L0_D, beta_L0_L, beta_A0_L, beta_L_D, beta_A0_D,
                        eta= rep(0.1,3), nu = rep(1.1,3), tau = 1){

  Delta <- ID <- Time_Disease <- Time <- NULL

  # Matrix for results
  res <- matrix(nrow = B, ncol = 8)
  colnames(res) <- c("Est 0", "Est 1", "SE 0", "SE 1", "Lower 0", "Lower 1",
                     "Upper 0", "Upper 1")

  for(b in 1:B){
    # Generate data
    data_boot <- simDisease(N = N, eta = eta, nu = nu, beta_L0_D = beta_L0_D, beta_L0_L = beta_L0_L,
                        beta_A0_L = beta_A0_L, beta_L_D = beta_L_D, beta_A0_D = beta_A0_D, cens = 1)

    # Disease events
    Disease_events <- data_boot[Delta == 2]

    # Disease people
    Disease_peeps <- data_boot[ID %in% Disease_events$ID]

    # Setting T_0 to debut time of diabetes
    Disease_peeps[, Time_Disease := Time - min(Time), by = ID]

    # Removing the new Time 0
    Disease_peeps <- Disease_peeps[Delta != 3]

    # Kaplan meyer fit
    boot_fit <- survival::survfit(Surv(Time_Disease, Delta == 1) ~ A0, data = Disease_peeps)

    # Save estimate of survival probability, SE and confidence limits
    preds <- summary(boot_fit, times = tau, data.frame = TRUE)

    # The KM estimate, SE's, lower CI, and upper CI is returned
    res[b, ] <- (c(preds[,5], preds[,7], preds[,9], preds[,10]))
  }

  return(res)
}
