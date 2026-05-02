#' Niche Divergence Plane Analysis for Two Species
#'
#' @param data_a Presence-background data for Species A
#' @param data_b Presence-background data for Species B
#' @param formula Formula for the glmnet_mx models
#' @param variables Vector of variable names to analyze
#' @export
ndp_analysis <- function(data_a, data_b, formula, variables) {
  # Fit models
  model_a <- glmnet_mx(data_a$p, data_a[, -1], formula)
  model_b <- glmnet_mx(data_b$p, data_b[, -1], formula)

  results <- list()

  for (var in variables) {
    spa <- get_response_curve(model_a, var)
    spb <- get_response_curve(model_b, var)
    ov <- beta_overlap(spa, spb)

    diss <- niche_diss(spa, spb, ov)
    excl <- niche_excl(spa, spb)

    results[[var]] <- data.frame(variable = var, dissimilarity = diss, exclusivity = excl)
  }

  do.call(rbind, results)
}

#' Bootstrap Niche Divergence Plane
#'
#' @param data_a Presence-background data for Species A
#' @param data_b Presence-background data for Species B
#' @param formula Formula for the glmnet_mx models
#' @param variables Vector of variable names to analyze
#' @param n_iter Number of bootstrap iterations
#' @export
ndp_bootstrap <- function(data_a, data_b, formula, variables, n_iter = 100) {
  boot_results <- list()

  for (i in 1:n_iter) {
    # Sample presence and background separately to maintain proportions?
    # Or just sample with replacement from the whole dataset.

    idx_a <- sample(1:nrow(data_a), replace = TRUE)
    idx_b <- sample(1:nrow(data_b), replace = TRUE)

    res <- ndp_analysis(data_a[idx_a, ], data_b[idx_b, ], formula, variables)
    res$iteration <- i
    boot_results[[i]] <- res
  }

  do.call(rbind, boot_results)
}

#' Permutation Test for Niche Divergence Plane
#'
#' @param data_a Presence-background data for Species A
#' @param data_b Presence-background data for Species B
#' @param formula Formula for the glmnet_mx models
#' @param variables Vector of variable names to analyze
#' @param n_iter Number of permutation iterations
#' @export
ndp_permutation <- function(data_a, data_b, formula, variables, n_iter = 100) {
  combined_data <- rbind(data_a, data_b)
  n_a <- nrow(data_a)
  n_total <- nrow(combined_data)

  perm_results <- list()

  for (i in 1:n_iter) {
    shuffled_idx <- sample(1:n_total)
    perm_a <- combined_data[shuffled_idx[1:n_a], ]
    perm_b <- combined_data[shuffled_idx[(n_a + 1):n_total], ]

    res <- ndp_analysis(perm_a, perm_b, formula, variables)
    res$iteration <- i
    perm_results[[i]] <- res
  }

  do.call(rbind, perm_results)
}
