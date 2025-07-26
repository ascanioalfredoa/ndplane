# Project Overview

- Project Name: ndplane
- Description: An R package to compare ecological niches using the Niche Divergence Plane methodology and framework. It enables quantification and classification of niche divergence scenarios across multiple niche axes.
- Key Technologies: R (>= 4.1.0), dplyr, furrr, future, pracma, progressr, readr, rlang, terra, tibble, tidyr, virtualspecies, ggplot2, purrr, progress, sf
- Architecture: Standard R package structure with source code in R/ directory, tests in tests/, and additional domain-specific analysis code in subfolders such as NDP_Salamanders, NDP_Simulations, NDP_VirtualSpecies.

# Getting Started

## Prerequisites
- R version 4.1.0 or higher
- Recommended: R development environment (RStudio)
- devtools package for development and installation

## Installation
```r
# Install development version from GitHub
install.packages("devtools") # if devtools not installed
devtools::install_github("ascanioalfredoa/ndplane")
```

## Basic Usage
```r
library(ndplane)
# Use ndplane functions to perform niche divergence analyses
```

## Running Tests
- Run tests using testthat framework:
```r
library(testthat)
test_package("ndplane")
```

# Project Structure

- R/: Core R source code files for the ndplane package
- tests/: Unit tests implemented with testthat
- README.Rmd, README.md: Project overview and basic usage
- DESCRIPTION: Package metadata and dependencies
- man/: Documentation files (help pages)
- NDP_Salamanders/: Scripts and data for Salamanders niche divergence analyses
- NDP_Simulations/: Simulation scripts and results
- NDP_VirtualSpecies/: Supporting scripts, data, and batch scripts for virtual species analyses

# Development Workflow

- Follow standard R package development practices
- Use roxygen2 comments in R/ source files to document functions
- Test with testthat framework
- Build and check package with devtools devtools::check()
- Use README.Rmd for dynamic README generation

# Key Concepts

- Niche Divergence Plane methodology: framework for ecological niche comparison
- Ecological niche axes: multiple environmental variables or niche dimensions analyzed
- Supports simulation and empirical data analyses

# Common Tasks

- Install package dependencies: use install.packages() as needed
- Run simulations or analyses: execute R scripts in project subdirectories
- Develop new features: add R scripts to R/ and document with roxygen2
- Test code thoroughly using tests/testthat/

# Troubleshooting

- Check R version compatibility
- Ensure all dependencies are installed
- Use devtools::check() to diagnose problems
- Run tests to identify failing code

# References

- CRAN R packages used: https://cran.r-project.org/
- devtools: https://cran.r-project.org/package=devtools
- testthat: https://testthat.r-lib.org/
- roxygen2 documentation: https://roxygen2.r-lib.org/

---

*Note: Some assumptions were made based on typical R package structures and content listed. Please verify details especially concerning project-specific workflows and domain concepts.*
