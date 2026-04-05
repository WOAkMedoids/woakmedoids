# WOAkMedoids 0.2.2

## Added
- Reimplemented performance-critical routines in C++ via Rcpp.
- Parallel fitness evaluation using RcppParallel.
- Early stopping mechanism to accelerate convergence.

## Changed
- Fused cluster assignment and distance accumulation into a single loop for improved performance.
- Refactored internal WOA implementation for better maintainability.

## Fixed
- Minor internal fixes.

# WOAkMedoids 0.1.0

## Added
- Initial pure-R implementation of WOA-kMedoids.
