# Handle global variables for R CMD check
utils::globalVariables(c("a", "b", ":="))

#' Normalize Suitability Raster from Virtual Species
#'
#' Normalizes a habitat suitability raster from a `virtualspecies` object so that all pixel values sum to 1.
#'
#' @param vsp A `virtualspecies` object generated with `generateSpFromFun()`.
#'
#' @return A `SpatRaster` object with normalized values.
#' @export
normalize_suitability <- function(vsp) {
    vsp_ras <- vsp$suitab.raster
    norm_vsp <- vsp_ras / sum(terra::values(vsp_ras), na.rm = TRUE)
    return(norm_vsp)
}

#' Calculate Niche Similarity Metrics (Schoener's D and Hellinger's I)
#'
#' Computes niche similarity metrics based on two normalized suitability rasters.
#'
#' @param sp1 A normalized suitability `SpatRaster` (species A).
#' @param sp2 A normalized suitability `SpatRaster` (species B).
#' @param method Character vector specifying one or both of: "D" (Schoener's D) and "I" (Hellinger's I).
#'
#' @return A data frame with columns for D and I (NA if not computed).
#' @export
calculate_similarity <- function(sp1, sp2, method = c("D", "I")) {
    method <- match.arg(method, several.ok = TRUE)

    result <- list(D = NA, I = NA)

    if ("D" %in% method) {
        result$D <- 1 - 0.5 * sum(abs(terra::values(sp1) - terra::values(sp2)), na.rm = TRUE)
    }
    if ("I" %in% method) {
        result$I <- 1 - 0.5 * sqrt(sum((sqrt(terra::values(sp1)) - sqrt(terra::values(sp2)))^2, na.rm = TRUE))
    }

    return(as.data.frame(result))
}

#' Define Variable-Specific Parameter Spaces for Virtual Species Simulations
#'
#' Builds a list of parameter grids (a, b, alpha, gamma) for each environmental variable used in beta response simulations. Optionally generates a matching raster stack with values sampled within each variable's domain.
#'
#' @param param_ranges A named list where each element is a list of vectors: `a`, `b`, `alpha`, `gamma` for a variable.
#' @param raster_stack Optional `SpatRaster` with layers named as in `param_ranges`. If NULL, random rasters are created.
#' @param raster_dim Dimensions (rows, cols) of random rasters if `raster_stack` is not provided (default = c(30, 30)).
#' @param save_to Optional folder path to save each variable's parameter grid as a .csv file.
#'
#' @return A list with:
#' \describe{
#'   \item{env_raster}{A `SpatRaster` stack of environmental layers}
#'   \item{parameter_grids}{A named list of parameter `tibble`s per variable}
#' }
#' @export
define_vsp_param_space <- function(
        param_ranges,
        raster_stack = NULL,
        raster_dim = c(30, 30),
        save_to = NULL
) {

    var_names <- names(param_ranges)

    # 1. Build or validate raster stack
    if (is.null(raster_stack)) {
        nrows <- raster_dim[1]
        ncols <- raster_dim[2]
        env_stack <- NULL

        for (var in var_names) {
            r_params <- param_ranges[[var]]
            a_range <- min(r_params$a)
            b_range <- max(r_params$b)

            r <- terra::rast(nrows = nrows, ncols = ncols,
                             vals = stats::runif(nrows * ncols, min = a_range, max = b_range))
            names(r) <- var

            if (is.null(env_stack)) {
                env_stack <- r
            } else {
                env_stack <- c(env_stack, r)
            }
        }
    } else {
        env_stack <- raster_stack
        if (!all(var_names %in% names(env_stack))) {
            stop("Raster stack is missing layers for: ", paste(setdiff(var_names, names(env_stack)), collapse = ", "))
        }
    }

    # 2. Generate parameter grids
    param_grids <- list()
    for (var in var_names) {
        r <- param_ranges[[var]]
        grid <- expand.grid(a = r$a, b = r$b, alpha = r$alpha, gamma = r$gamma)
        grid <- dplyr::filter(grid, a < b)
        grid <- tibble::as_tibble(grid)
        param_grids[[var]] <- grid

        if (!is.null(save_to)) {
            if (!dir.exists(save_to)) dir.create(save_to, recursive = TRUE)
            readr::write_csv(grid, file = file.path(save_to, paste0("params_", var, ".csv")))
        }
    }

    # 3. Output
    return(list(
        env_raster = env_stack,
        parameter_grids = param_grids
    ))
}

#' Randomly sample one parameter set per variable
#'
#' @param param_grids Named list of parameter data frames (from define_vsp_param_space)
#' @param strategy Sampling strategy: "random" or "first"
#'
#' @return Named list of sampled rows as named vectors
#' @export
sample_params <- function(param_grids, strategy = "random") {
    lapply(param_grids, function(grid) {
        if (strategy == "random") {
            grid[sample(nrow(grid), 1), ]
        } else {
            grid[1, ]
        }
    })
}

#' Rescaled beta response function (from virtualspecies internal)
#'
#' This function returns habitat suitability values based on a rescaled beta function.
#' It is used in virtual species simulations to define response curves along environmental gradients.
#'
#' @param x Numeric. Vector of environmental values.
#' @param p1 Numeric. Lower limit of species tolerance.
#' @param p2 Numeric. Upper limit of species tolerance.
#' @param alpha Numeric. Shape parameter for the rising edge.
#' @param gamma Numeric. Shape parameter for the falling edge.
#'
#' @return A numeric vector of habitat suitability values (0 to max k).
#' @export
vsp_betaFun <- function(x, p1, p2, alpha, gamma) {
    k <- 1 / ((alpha * (p2 - p1) / (alpha + gamma))^alpha) /
        ((gamma * (p2 - p1) / (alpha + gamma))^gamma)

    ifelse(
        x > p1 & x < p2,
        k * ((x - p1)^alpha) * ((p2 - x)^gamma),
        0
    )
}


#' Format parameter list for virtualspecies::formatFunctions
#'
#' @param param_list A named list of parameter rows (each 1-row tibble)
#'
#' @return A list formatted for generateSpFromFun()
#' @export
format_vsp <- function(param_list) {
    args <- purrr::map(param_list, function(p, var) {
        stats::setNames(
            c("vsp_betaFun", p$a, p$b, p$alpha, p$gamma),
            c("fun", "p1", "p2", "alpha", "gamma"))
    }, var = names(param_list))

    names(args) <- names(param_list)  # preserve variable names for formatFunctions()

    rlang::eval_bare(rlang::call2(virtualspecies::formatFunctions, !!!args))

}

#' Calculate niche divergence indices from betaPDF
#'
#' Computes the components of the Niche Divergence Plane (NDP):
#' - ND: Niche Dissimilarity (Y-axis)
#' - NE: Niche Exclusivity (X-axis)
#' - ND_magnitude: Euclidean magnitude in NDP space
#' - ND_angle: Direction of divergence in degrees (atan2(ND, NE))
#'
#' @param param_A Named vector with a, b, alpha, gamma
#' @param param_B Named vector with a, b, alpha, gamma
#'
#' @return A tibble with ND, NE, ND_magnitude, ND_angle
#' @export
nd_indices <- function(param_A, param_B) {
    A <- betaPDF(a = param_A$a, b = param_A$b, alpha = param_A$alpha, gamma = param_A$gamma)
    B <- betaPDF(a = param_B$a, b = param_B$b, alpha = param_B$alpha, gamma = param_B$gamma)
    ov <- beta_overlap(A, B)

    ND <- round(niche_diss(A, B, ov), 3)
    NE <- round(niche_excl(A, B), 3)
    ND_magnitude <- sqrt(ND^2 + NE^2)
    ND_angle <- atan2(ND, NE) * (180 / pi)  # ND is Y, NE is X

    tibble::tibble(
        ND = ND,
        NE = NE,
        ND_magnitude = ND_magnitude,
        ND_angle = ND_angle
    )
}
#' Simulate two virtual species and compute similarity metrics
#'
#' @param params_A Formatted functions list for species A
#' @param params_B Formatted functions list for species B
#' @param env_stack SpatRaster stack used to generate species
#'
#' @return List with rasters and similarity indices
#' @export
simulate_vsp <- function(params_A, params_B, env_stack) {
    vsp_A <- virtualspecies::generateSpFromFun(
        raster.stack = env_stack,
        parameters = params_A,
        formula = paste(names(params_A), collapse = " + "),
        plot = FALSE
    )
    vsp_B <- virtualspecies::generateSpFromFun(
        raster.stack = env_stack,
        parameters = params_B,
        formula = paste(names(params_B), collapse = " + "),
        plot = FALSE
    )

    sp1 <- normalize_suitability(vsp_A)
    sp2 <- normalize_suitability(vsp_B)

    sim <- calculate_similarity(sp1, sp2, method = c("D", "I"))

    list(D = sim$D, I = sim$I, sp1 = sp1, sp2 = sp2)
}

#' Run NDP simulations using virtual species and parameter space
#'
#' @param vsp_space Output of define_vsp_param_space()
#' @param n_comparisons Number of pairwise comparisons to simulate
#' @param strategy Sampling strategy: "random" or "first"
#' @param seed Optional seed for reproducibility
#' @param progress Logical; if TRUE, show progress bar
#' @param save_as Optional file path (.csv or .rds) to save the results
#'
#' @return A tibble with NDP and similarity metrics for each variable and comparison
#' @export
run_vsp_sim <- function(vsp_space, n_comparisons = 10, strategy = "random",
                        seed = NULL, progress = TRUE, save_as = NULL) {
    if (!requireNamespace("progress", quietly = TRUE) && progress) {
        message("progress package not installed; disabling progress bar.")
        progress <- FALSE
    }

    if (!is.null(seed)) set.seed(seed)

    param_grids <- vsp_space$parameter_grids
    env_stack <- vsp_space$env_raster
    vars <- names(param_grids)

    if (progress) pb <- progress::progress_bar$new(
        format = "Simulating [:bar] :percent eta: :eta",
        total = n_comparisons, clear = FALSE, width = 60
    )

    results <- purrr::map_dfr(1:n_comparisons, function(i) {
        if (progress) pb$tick()

        params_A_raw <- sample_params(param_grids, strategy)
        params_B_raw <- sample_params(param_grids, strategy)

        formatted_A <- format_vsp(params_A_raw)
        formatted_B <- format_vsp(params_B_raw)

        sim <- simulate_vsp(formatted_A, formatted_B, env_stack)

        purrr::map_dfr(vars, function(var) {
            nd <- nd_indices(params_A_raw[[var]], params_B_raw[[var]])
            dplyr::mutate(
                nd,
                comparison_id = i,
                variable = var,
                D = sim$D,
                I = sim$I,
                a_A = params_A_raw[[var]]$a,
                b_A = params_A_raw[[var]]$b,
                alpha_A = params_A_raw[[var]]$alpha,
                gamma_A = params_A_raw[[var]]$gamma,
                a_B = params_B_raw[[var]]$a,
                b_B = params_B_raw[[var]]$b,
                alpha_B = params_B_raw[[var]]$alpha,
                gamma_B = params_B_raw[[var]]$gamma
            )
        })
    })

    if (!is.null(save_as)) {
        if (grepl("\\.csv$", save_as)) {
            readr::write_csv(results, save_as)
        } else if (grepl("\\.rds$", save_as)) {
            saveRDS(results, save_as)
        } else {
            warning("Unsupported file format in `save_as`. Must be .csv or .rds")
        }
    }

    return(results)
}

#' Plot results from run_vsp_sim
#'
#' Generate a customizable scatterplot from NDP simulation results (e.g., ND vs D).
#'
#' @param results A data frame or tibble returned by `run_vsp_sim()`
#' @param x Character. Column name for the x-axis (default = "ND")
#' @param y Character. Column name for the y-axis (default = "D")
#' @param color_by Optional character. Variable to use for color (default = "variable")
#' @param facet_by Optional character. Variable to use for faceting (default = NULL)
#' @param legend_title Optional string. Legend title to display (default = "variable")
#' @param mapped_args Optional list of mapped aesthetics (e.g., list(shape = "variable"))
#' @param fixed_args Optional list of fixed aesthetics (e.g., list(size = 2, alpha = 0.6))
#'
#' @return A `ggplot2` object
#' @export
plot_ndp_results <- function(results,
                             x = "ND",
                             y = "D",
                             color_by = "variable",
                             facet_by = NULL,
                             legend_title = "variable",
                             mapped_args = list(),
                             fixed_args = list()) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The 'ggplot2' package is required for this function.")
    }
    if (!requireNamespace("rlang", quietly = TRUE)) {
        stop("The 'rlang' package is required for tidy evaluation.")
    }

    results <- as.data.frame(results)

    # Build aes expression dynamically
    aes_expr <- rlang::expr(ggplot2::aes(
        x = !!rlang::sym(x),
        y = !!rlang::sym(y)
    ))

    if (!is.null(color_by)) {
        aes_expr <- rlang::call_modify(aes_expr, colour = rlang::sym(color_by))
    }

    for (arg in names(mapped_args)) {
        aes_expr <- rlang::call_modify(aes_expr, !!arg := rlang::sym(mapped_args[[arg]]))
    }

    # Combine into argument list for geom_point
    geom_args <- c(
        list(mapping = eval(aes_expr), inherit.aes = FALSE),
        fixed_args
    )

    p <- ggplot2::ggplot(results) +
        do.call(ggplot2::geom_point, geom_args)

    if (!is.null(facet_by)) {
        p <- p + ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet_by)))
    }

    if (!is.null(legend_title) && !is.null(color_by)) {
        p <- p + ggplot2::labs(color = legend_title)
    } else if (is.null(legend_title)) {
        p <- p + ggplot2::guides(color = "none")
    }

    p <- p + ggplot2::theme_minimal() +
        ggplot2::labs(x = x, y = y)

    return(p)
}
