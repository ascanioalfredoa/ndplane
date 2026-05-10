# Development Plan for ndplane R Package

This document outlines the roadmap for the development and release of the `ndplane` package, focusing on a lightweight, functional first version (v0.1.0) and the path towards a CRAN-ready release.

## Vision
To provide a focused, lightweight R package for quantifying and classifying ecological niche divergence using the Niche Divergence Plane (NDP) framework (Ascanio et al. 2024).

## Phase 1: Core Functionality (v0.1.0) - [COMPLETED]
- [x] **Lightweight Modeling**: Implement `glmnet_mx` to allow Maxent-like modeling without heavy dependencies (`maxnet`, `dismo`).
- [x] **NDP Indices**: Implement Niche Dissimilarity (`niche_diss`) and Niche Exclusivity (`niche_excl`) using Base R logic and a custom trapezoidal integration (`trapz`).
- [x] **Base R Plotting**: Provide high-quality visualizations of response curves and the NDP using Base R `graphics`.
- [x] **Uncertainty & Testing**: Implement `ndp_bootstrap` for uncertainty estimation and `ndp_permutation` for null hypothesis testing.
- [x] **Documentation**: Set up initial vignettes and a `pkgdown` website.

## Phase 2: Refinement and Cleaning
- [ ] **Code Tidying**: Ensure all internal code follows the Tidyverse style guide while remaining dependency-free where possible.
- [ ] **Argument Validation**: Add robust input checks to all exported functions to provide helpful error messages.
- [ ] **Performance Optimization**: Optimize the bootstrap and permutation loops for speed, potentially exploring Base R parallelization (`parallel` package).

## Phase 3: Package Ecosystem
- [ ] **Helper Packages**: Move non-core functions (spatial extraction, virtual species simulations, specialized plotting) to helper packages (e.g., `ndplane.spatial`, `ndplane.sim`).
- [ ] **Compatibility**: Ensure `ndplane` can easily ingest data from common ENM workflows (e.g., `terra`, `sf`, `ENMTools`).

## Phase 4: CRAN Readiness
- [ ] **Testing**: Achieve high test coverage for core logic and edge cases.
- [ ] **Manuals**: Complete all `.Rd` files with detailed documentation and examples.
- [ ] **CI/CD**: Ensure the package passes `R CMD check` on all major platforms (Linux, Windows, macOS).
- [ ] **Submission**: Submit to CRAN.

## Contribution Guidelines
- Favor Base R for core logic to keep dependencies to a minimum.
- All new features should be accompanied by unit tests.
- Documentation should be updated with every significant change.
