# cowfit: LMMs and GLMMs in genetic evaluation <img src="fig/cowfit.png" width="120" align="right" />

This package extends the functionality of the package `pedigreemm` and contains the following improvements

1. Resolves a critical bug in the prediction of random effects.
2. Allows for animal models with single observations.
3. Allows for random regression models.
4. Allows to prespecify the variance components for faster estimation of random effects.

Also, the same models can be fitted in the **Bayesian framework**. The implementation is based on the package `brms` with some adaptations

1. Uses the same syntax for defining correlation between animals as `pedigreemm` (i.e. with argument `pedigree = list(...)`).
2. Avoids the calculation of the nummerator relationship matrix A (Covariance matrix between random effects).
3. Allows to prespecify the variance components for faster estimation of random effects.

## Installation

Install the package `devtools` (with `install.packages("devtools")`) and run the following line of code
```
devtools::install_github("retodomax/cowfit", dependencies=TRUE)
```

## Documentation

* [Getting started with cowfit](https://retodomax.github.io/cowfit/getting_started)
* [Theory]
