
<!-- README.md is generated from README.Rmd. Please edit that file -->

# simevent

<!-- badges: start -->
<!-- badges: end -->

The goal of simevent is to provide functions for the generation and
analysis of complex continuous time health care data.

## Installation

You can install the development version of simevent from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("miclukacova/simevent")
```

## Example 1: simEventData

The underlying function is called `simEventData`. One can specify
various arguments

``` r
library(simevent)
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
#>       ID       Time Delta       L0     L    A0     A
#>    <int>      <num> <num>    <num> <num> <int> <num>
#> 1:     1 0.09302029     3 69.34624     0     1     0
#> 2:     1 1.15310728     1 69.34624     1     1     0
#> 3:     2 0.68901160     2 69.89826     0     0     0
#> 4:     2 2.72769569     0 69.89826     0     0     1
#> 5:     3 1.78458962     1 37.39783     0     1     0
#> 6:     4 0.22118200     2 32.73879     0     0     0
```

One can visualize the data by

``` r
plotEventData(data)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

## Example 2: Survival Data

## Example 3: Competing Risk Data
