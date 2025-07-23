# birdie 0.7.0

* Support for 2020 decennial census data (#19)
* Fix bug preventing `p_r="estimate"` in `bisg_me()` (#18)
* Update CITATION
* Add `rstan` to `Suggests` due to its use in `Makefile` (#21)

# birdie 0.6.1

* Switch to an inverse-gamma prior for the random effects scale parameter in 
  the mixed model.
* Add Jacobian adjustment to mixed model M step to improve accuracy and reduce
  prior sensitivity.

# birdie 0.5.0

* Gibbs sampler for Categorical-Dirichlet and Normal linear models
* `simulate.bisg()` to perform multiple imputation from BISG probabilities

# birdie 0.4.0

* Support for continuous outcomes
* Support for observation/likelihood weights
* New method for specifying outcome model type
* New default empirical Bayes prior for Categorical-Dirichlet model

# birdie 0.3.0

* New interfaces and implementations for BISG
* New computational backend and interface for BIRDiE models
* Updated documentation and examples
* Unit tests for accuracy and usability

# birdie 0.1.0

* Basic BISG functionality and initial model support
