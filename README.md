
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simevent

<!-- badges: start -->
<!-- badges: end -->

The goal of `simevent` is to provide functions for the generation and
analysis of complex continuous time health care data.

## Installation

You can install the development version of simevent from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("miclukacova/simevent")
```

## Example 1: simEventData

We load the package

``` r
library(simevent)
```

The underlying function is called `simEventData`. One can specify
various arguments, as for example the `N` and `beta` argument. The
`N`argument lets the user specify number of individuals in the
simulation. The `beta` argument lets the user specify the effects of
processes and covariates on the intensities of the processes. We define
the various arguments

``` r
# Number of individuals
N <- 100
# Effect on event 0 
beta0 <- c(0, 0, 0, 0)
# Effect on event 1
beta1 <- c(1, -1, 1, -1)
# Effect on event 2 (A)
beta2 <- c(0, -1, 0, 0.5)
# Effect on event 3 (L)
beta3 <- c(0, 0, 1, 0)
beta <- cbind(beta0, beta1, beta2, beta3)
```

And then call the function

``` r
data <- simEventData(N = N, beta = beta)
```

The simulated data looks like

``` r
head(data)
#> Key: <ID>
#>       ID      Time Delta       L0     L    A0     A
#>    <int>     <num> <num>    <num> <num> <int> <num>
#> 1:     1  8.400588     2 56.93411     0     1     0
#> 2:     1 10.302457     0 56.93411     0     1     1
#> 3:     2  1.845571     2 43.15434     0     0     0
#> 4:     2  6.482616     1 43.15434     0     0     1
#> 5:     3  5.606393     0 51.53736     0     1     0
#> 6:     4  2.760928     1 41.77216     0     0     0
```

One can visualize the data by

``` r
plotEventData(data)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

### Classical Estimation

We simulate for $N = 10^4$ individuals

``` r
data <- simEventData(N = 10^4, beta = beta)
```

In order to fit a Cox proportional hazards model with the `survival`
package, the data needs to be transformed into the so called format,
this can be done by the function \`IntFormatData.\`\`\` Furthermore two
at risk indicators need to be created, one that indicates when the
individual is at risk for event $A$ and one that indicates that the
individual is at risk for event $L$.

``` r
data_int <- IntFormatData(data)
data_int[, at_risk_2 := as.numeric(A == 0)]
data_int[, at_risk_3 := as.numeric(L == 0)]
```

Data in the format looks like

``` r
head(data_int)
#> Key: <ID>
#>       ID     Time Delta       L0     L    A0     A     k   tstart    tstop
#>    <int>    <num> <num>    <num> <num> <int> <num> <int>    <num>    <num>
#> 1:     1 1.022917     3 60.88006     0     1     0     1 0.000000 1.022917
#> 2:     1 1.097404     1 60.88006     1     1     0     2 1.022917 1.097404
#> 3:     2 3.025314     3 35.68031     0     1     0     1 0.000000 3.025314
#> 4:     2 6.343352     0 35.68031     1     1     0     2 3.025314 6.343352
#> 5:     3 0.900248     0 44.08172     0     1     0     1 0.000000 0.900248
#> 6:     4 1.596860     0 37.91036     0     1     0     1 0.000000 1.596860
#>    at_risk_2 at_risk_3
#>        <num>     <num>
#> 1:         1         1
#> 2:         1         0
#> 3:         1         1
#> 4:         1         0
#> 5:         1         1
#> 6:         1         1
```

Cox proportional hazards models for the death process and operation
process can be fitted by the following code

``` r
library(survival)
# Death process
survfit_death <- coxph(Surv(tstart, tstop, Delta == 1) ~ I(L0/50) + A0 + L + A, 
                       data = data_int)
# Operation process
survfit_oper <- coxph(Surv(tstart, tstop, Delta == 2) ~ I(L0/50) + A0 + L, 
                      data = data_int[at_risk_2 == 1])
```

The regression results can be seen by the summary call

``` r
survfit_death |> summary()
#> Call:
#> coxph(formula = Surv(tstart, tstop, Delta == 1) ~ I(L0/50) + 
#>     A0 + L + A, data = data_int)
#> 
#>   n= 14443, number of events= 6286 
#> 
#>              coef exp(coef) se(coef)      z Pr(>|z|)    
#> I(L0/50)  0.94596   2.57529  0.05484  17.25   <2e-16 ***
#> A0       -0.99750   0.36880  0.02716 -36.73   <2e-16 ***
#> L         0.97131   2.64139  0.02972  32.69   <2e-16 ***
#> A        -0.95891   0.38331  0.04183 -22.92   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>          exp(coef) exp(-coef) lower .95 upper .95
#> I(L0/50)    2.5753     0.3883    2.3128    2.8675
#> A0          0.3688     2.7115    0.3497    0.3890
#> L           2.6414     0.3786    2.4919    2.7998
#> A           0.3833     2.6089    0.3531    0.4161
#> 
#> Concordance= 0.673  (se = 0.004 )
#> Likelihood ratio test= 2645  on 4 df,   p=<2e-16
#> Wald test            = 2640  on 4 df,   p=<2e-16
#> Score (logrank) test = 2733  on 4 df,   p=<2e-16
```

## Example 2: Survival Data

You can simulate data from a survival setting with the function
`simSurvData`.

``` r
data <- simSurvData(100)
plotEventData(data, title = "Survival Data")
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

One can again specify the effects of $A_0$ and $L_0$ on the risk of
death and censoring by the `beta` argument.

``` r
# No effect of L0 and A0 on censoring process
beta_C <- c(0,0)
# Effect of L0 and A0 on death process
beta_D <- c(1,-1)

beta <- cbind(beta_C, beta_D)
```

And specify the parameters of the Weibull intensity for the censoring
and death process.

``` r
eta <- c(0.2, 0.2)
nu <- c(1.05, 1.05)
```

We now call the function and visualize the data

``` r
data <- simSurvData(100, beta = beta, eta = eta, nu = nu)
plotEventData(data, title = "Survival Data")
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

## Example 3: Competing Risk Data

You can simulate data from a competing risk setting with the function
`simCRdata`. The arguments `beta`, `eta`, `nu`, work in a similar maner
as above.

``` r
data <- simCRdata(100)
plotEventData(data, title = "Competing Risk Data")
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

## Example 4: Type 2 Diabetes

The function `simT2D` simulates health care data from a setting where
patients can experience $3$ different events: Censoring (0), Death (1)
and Type-2-Diabetes (2). The various arguments allow for the different
scenarios, and you can read about them on the help page

``` r
?simT2D
```

Below is a function call to `simT2D`

``` r
data <- simT2D(N = 100,
               sex = FALSE, 
               cens = 1,
               eta = c(0.1,0.3,0.1,0.1), 
               nu = c(1.1,1.3,1.1,1.1),
               beta_L0_L = 1, 
               beta_A0_L = -1.1, 
               beta_L_D = 1, 
               beta_L0_D = 0)

plotEventData(data, title = "T2D data")
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="100%" />

## Example 5: Unobserved Covariate Setting

The `simConfounding` function was created to simulate data from a
setting where we have an unobserved confounding process. You can read
about the function on the help page

``` r
?simConfounding
```

One can simulate data from the default setting by the function call

``` r
data <- simConfounding(100)
```

And one can simulate from user specified scenarios by the function call

``` r
data <- simConfounding(N = 100,
                       beta_L_A = 1,
                       beta_L_D = 1,
                       beta_A_D = -1,
                       beta_A_L = -0.5,
                       beta_L0_A = 1,
                       eta = rep(0.1, 4),
                       nu = rep(1.1, 4),
                       followup = 5,
                       cens = 1,
                       op = 1)
```

For example the function call above simulates from a setting with the
operation/treatment event (op = 1), where there is a censoring process
(cens = 1), and where after 5 time units everybody is censored (followup
= 5). We can again visulize data by the function call to
`plotEventData`.

``` r
plotEventData(data, title = "Confounding setting")
```

<img src="man/figures/README-unnamed-chunk-20-1.png" width="100%" />
