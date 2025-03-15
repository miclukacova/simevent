#' Function that performs a simulation study of the Kaplan-Meier estimator on data
#' from a T2D diabetes setting. Data is simulated by use of the function simT2D.
#' Right censoring is incorporated, and the L1 variable is omitted. The estimand is
#' \deqn{S(\tau | A_0 = 1) = P(T \geq \tau | A_0 = 1) , \quad S(\tau | A_0 = 0) = P(T \geq \tau | A_0 = 0)}
#' For some \eqn{\tau \in \mathbb{R}}. \eqn{T} is here the survival time *post* T2D diagnose.
#' The survival function will be estimated by the Kaplan-Meier estimator.
#' Post T2D diagnose, the setting is simpler, as there are no competing events, only death and censoring.
#' The Kaplan-Meier estimator will then be an unbiased estimator of the survival probability.
#' The function generates B data sets with N individuals. For each data set we find the
#' Kaplan-Meier estimate of \eqn{P(T \geq \tau | A_0 = a_0)}, the SE, and the upper and
# lower confidence bands, for \deqn{a_0 \in \{0,1\}}.
#' @title simStudyT2D
#'
#' @param N A double of the number of individuals
#' @param B A double for the number of repetitions
#' @param beta_L0_D Specifies how L0 affects risk of death.
#' @param beta_L0_L Specifies how L0 affects risk of T2D.
#' @param beta_A0_L Specifies how A0 affects risk of T2D.
#' @param beta_L_D Specifies how change in the covariate process affects risk of death.
#' @param beta_A0_D Specifies how A0 affects risk of death.
#' @param eta Vector of length 4 of shape parameters for the Weibull intensity with parameterization
#' \deqn{\eta \nu t^{\nu - 1}}. Default is set to 0.1 for all events.
#' @param nu Vector of scale parameters for the Weibull hazard. Default is set to 1.1 for all events.
#' @param tau The time at which the survival probability is estimated.
#'
#' @returns Results of simulation study in the form of a B times 8 matrix. The columns are the estimates of
#' survival probability, SE, lower and upper CI for respectively A0 = 0 and A0 = 1.
#' @export
#'
#' @examples
#' simStudyT2D(N = 100, B = 10, beta_L0_D = 1, beta_L0_L = 1, beta_A0_L = -1,
#'             beta_L_D = 1, beta_A0_D = 0, tau = 1)
simStudyT2D <- function(N, B, beta_L0_D, beta_L0_L, beta_A0_L, beta_L_D, beta_A0_D,
                        eta= rep(0.1,4), nu = rep(1.1,4), tau = 1){

  Delta <- ID <- Time_T2D <- Time <- NULL

  # Matrix for results
  res <- matrix(nrow = B, ncol = 8)
  colnames(res) <- c("Est 0", "Est 1", "SE 0", "SE 1", "Lower 0", "Lower 1",
                     "Upper 0", "Upper 1")

  for(b in 1:B){
    # Generate data
    data_boot <- simT2D(N = N, eta = eta, nu = nu, sex = FALSE, beta_L0_D = beta_L0_D, beta_L0_L = beta_L0_L,
                        beta_A0_L = beta_A0_L, beta_L_D = beta_L_D, beta_A0_D = beta_A0_D, cens = 1)

    # T2D events
    T2D_events <- data_boot[Delta == 3]

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
