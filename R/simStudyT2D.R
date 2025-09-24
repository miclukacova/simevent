#' Simulation Study of the Kaplan-Meier Estimator Post T2D Diagnosis
#'
#' This function performs a simulation study evaluating the Kaplan-Meier estimator
#' in a Type 2 Diabetes (T2D) context. Data is generated using `simT2D`, with
#' right censoring incorporated. The covariate `L1` is omitted. The estimand is:
#'
#' \deqn{S(\tau | A_0 = a_0) = P(T \geq \tau | A_0 = a_0), \quad a_0 \in \{0, 1\}, \ \tau \in \mathbb{R}}
#'
#' Here, \eqn{T} is the survival time post T2D diagnosis. After T2D onset, only death and censoring occur (no competing risks),
#' making the Kaplan-Meier estimator unbiased for the survival probability.
#'
#' For each of the \code{B} simulation runs, a dataset with \code{N} individuals is generated. The Kaplan-Meier estimate,
#' standard error (SE), and 95% confidence interval (CI) for the survival function at time \code{tau} are computed
#' separately for \eqn{A_0 = 0} and \eqn{A_0 = 1}.
#'
#' @title simStudyT2D
#'
#' @param N Integer. Number of individuals in each simulated dataset.
#' @param B Integer. Number of simulation repetitions.
#' @param beta_L0_D Numeric. Effect of L0 on risk of death.
#' @param beta_L0_L Numeric. Effect of L0 on risk of T2D.
#' @param beta_A0_L Numeric. Effect of A0 on risk of T2D.
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
#' simStudyT2D(N = 100, B = 10, beta_L0_D = 1, beta_L0_L = 1, beta_A0_L = -1,
#'             beta_L_D = 1, beta_A0_D = 0, tau = 1)
simStudyT2D <- function(N, B, beta_L0_D, beta_L0_L, beta_A0_L, beta_L_D, beta_A0_D,
                        eta= rep(0.1,3), nu = rep(1.1,3), tau = 1){

  Delta <- ID <- Time_T2D <- Time <- NULL

  # Matrix for results
  res <- matrix(nrow = B, ncol = 8)
  colnames(res) <- c("Est 0", "Est 1", "SE 0", "SE 1", "Lower 0", "Lower 1",
                     "Upper 0", "Upper 1")

  for(b in 1:B){
    # Generate data
    data_boot <- simT2D(N = N, eta = eta, nu = nu, beta_L0_D = beta_L0_D, beta_L0_L = beta_L0_L,
                        beta_A0_L = beta_A0_L, beta_L_D = beta_L_D, beta_A0_D = beta_A0_D, cens = 1)

    # T2D events
    T2D_events <- data_boot[Delta == 2]

    # T2D people
    T2D_peeps <- data_boot[ID %in% T2D_events$ID]

    # Setting T_0 to debut time of diabetes
    T2D_peeps[, Time_T2D := Time - min(Time), by = ID]

    # Removing the new Time 0
    T2D_peeps <- T2D_peeps[Delta != 3]

    # Kaplan meyer fit
    boot_fit <- survival::survfit(Surv(Time_T2D, Delta == 1) ~ A0, data = T2D_peeps)

    # Save estimate of survival probability, SE and confidence limits
    preds <- summary(boot_fit, times = tau, data.frame = TRUE)

    # The KM estimate, SE's, lower CI, and upper CI is returned
    res[b, ] <- (c(preds[,5], preds[,7], preds[,9], preds[,10]))
  }

  return(res)
}
