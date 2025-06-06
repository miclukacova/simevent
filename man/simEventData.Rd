% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simEventData.R
\name{simEventData}
\alias{simEventData}
\title{Simulate Event Data}
\usage{
simEventData(
  N,
  beta = NULL,
  eta = NULL,
  nu = NULL,
  at_risk = NULL,
  term_deltas = c(0, 1),
  max_cens = Inf,
  add_cov = NULL,
  override_beta = NULL,
  max_events = 10,
  lower = 10^(-15),
  upper = 200,
  gen_A0 = NULL
)
}
\arguments{
\item{N}{A double for the number of simulated individuals}

\item{beta}{A matrix of doubles for the effects of covariates and events on the
intensities. The columns represent the events N0, N1, .... In the default case 4 events.
The rows represent covariates and processes: the first two rows determine the
effects of the baseline covariate \eqn{L0} and \eqn{A0} on the processes. The
next rows determine the effects of the additional baseline covariates, by default
named L1, L2,... on the processes, followed by the effects of the processes N0,
N1,.... The \eqn{\beta} matrix is by  default set to 0.}

\item{eta}{Vector of shape parameters for the baseline intensity. Default is set
to 0.1 for all events.}

\item{nu}{Vector of scale parameters for the baseline intensity. Default is set
to 1.1 for all events.}

\item{at_risk}{At risk function. The \code{at_risk} function determines whether
an individual is at risk for a specific event. The function takes as input a vector
\code{events} of event counts. The function returns a vector of 0's and 1's indicating
which events the subject is at risk for. Default is set to a setting where you
are always at risk for all events.}

\item{term_deltas}{Terminal events. Default is set so that event 0 and 1 are
terminal events.}

\item{max_cens}{A maximum censoring time. By default set to infinity. If the event
time is larger than this maximal censoring time the event is set to 0 with prob
ability 1 and the event time is set to \code{max_cens}.}

\item{add_cov}{Named list of random generator functions for the distributions of
additional baseline covariates. The functions should take the number of observations
as input. By default set to NULL.}

\item{override_beta}{This argument is a named list. The argument has two applications.
One possibility is instead of specifying the whole beta matrix, the user can
specify the relevant entries, and the rest will by default be 0. Imagine you want
to specify the effect of L0 on N1 to be equal to 2, this could be done by
override_beta = list("L0" = c("N1" = 2))
In general if you want the effect \eqn{\beta_{x,y}= z}, you can specify the
override_beta argument as
override_beta = list("x" = c("y" = z))
Where x and y are names of processes or covariates.}

\item{max_events}{Number of maximal events per individual}

\item{lower}{Number of maximal events per individual}

\item{upper}{Number of maximal events per individual}

\item{gen_A0}{Function for generation of A0 covariate. Function of N (number
of individuals) and L0 (baseline covariate).}
}
\value{
data.table containing the simulated data. There is a column for ID, time
of event (Time), event (Delta), baseline covariate (L0), Baseline Treatment (A0),
the count of the various events: N1, N2, .... In case of additional covariates
these are included in the data as well, named L1, L2, ....
}
\description{
\code{simEventData} is a function to simulate continuous time to event data, e.g. observational
healthcare data. The number of events simulated corresponds to the length of the \eqn{\eta}
and \eqn{\nu} vector, as well as the number of columns in the \eqn{\beta} matrix. By
default 4 different types of events are simulated. The first two being terminal
processes. The simulations build upon a counting process framework, where the events
have intensities given by
\deqn{\lambda^x(t) = \lambda_0^x \exp(\beta^T_x L)}
L here represents covariates and event counts. The baseline intensity \eqn{\lambda_0^x}
is given by
\deqn{\lambda_0^x(t)=\eta^x \nu^x t^{\nu^x - 1}}
The intensities of the various events depend upon previous events and
baseline covariates through the \eqn{\beta} matrix. The pre specified \eqn{\nu}
and \eqn{\eta} parameters also influence the intensities of the various events
through the baseline intensity.
}
\examples{
simEventData(N = 10)
}
