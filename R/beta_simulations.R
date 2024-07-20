#' Niche Divergence Calculations for a given Parameter Space based on the beta function curve
#'
#' @param psp Parameter space created using the set_parameter_space() function
#'
#' @return A data frame containing niche divergence plane values of Dissimilarity and Exclusivity, along with the parameters of each species response
#' @export
#'
#' @examples
#' # Set parameter domains
#' a = b <- round(1:3, 3)
#' alpha = gamma <- round(1)
#'
#' # Set parameter space (combination of parameters for any beta function)
#' par_space <- set_parameter_space(a, b, alpha, gamma)
#'
#' #Run calc
#' calculate_nichediv_b1(psp = par_space)
#'
calculate_nichediv_b1 <- function(psp) {
    par_space <- psp
    row_combs <- combn(1:nrow(par_space), 2, simplify = F)

    results_names <- c("spA_par", "spB_par", "Dissimilarity", "Exclusivity",
                       paste(names(par_space), "_A", sep = ""),
                       paste(names(par_space), "_B", sep = "")
    )

    results <- dplyr::as_tibble(array(data = 0,
                                      dim = c(1, length(results_names)),
                                      dimnames = list(1, results_names)))

    res <- progressr::with_progress({
        p <- progressr::progressor(steps = length(row_combs))
        future::plan('future::sequential')
        row_combs |> furrr::future_map(\(x) {
            p()
            i <- x[1]; j <- x[2]
            spa <- betaPDF_ps(par_space[i, ])
            spb <- betaPDF_ps(par_space[j, ])
            ov <- beta_overlap(spa, spb)
            results[1, ] <- c(i, j,
                              niche_diss(spa, spb, ov),
                              niche_excl(spa, spb),
                              par_space[i, ], par_space[j, ])
            results
        }) |>
            dplyr::bind_rows()
    }
    )
    res
}


#' (Parallel) Niche Divergence Calculations for a given Parameter Space based on the beta function curve
#'
#' @param psp Parameter space created using the set_parameter_space() function
#' @param n number of cores to use in parallel processing (see future::plan())
#'
#' @return A data frame containing niche divergence plane values of Dissimilarity and Exclusivity, along with the parameters of each species response
#' @export
#'
#' @examples
#' # Set parameter domains
#' a = b <- round(1:3, 3)
#' alpha = gamma <- round(1)
#'
#' # Set parameter space (combination of parameters for any beta function)
#' par_space <- set_parameter_space(a, b, alpha, gamma)
#'
#' #Run calc
#' calculate_nichediv_bpar(psp = par_space)
#' #By default function uses 2 cores, but you can run
#' #n = parallelly::availableCores()-1
calculate_nichediv_bpar <- function(psp, n = 2) {
    par_space <- psp
    row_combs <- combn(1:nrow(par_space), 2, simplify = F)

    results_names <- c("spA_par", "spB_par", "Dissimilarity", "Exclusivity",
                       paste(names(par_space), "_A", sep = ""),
                       paste(names(par_space), "_B", sep = "")
    )

    results <- dplyr::as_tibble(array(data = 0,
                                      dim = c(1, length(results_names)),
                                      dimnames = list(1, results_names)))

    progressr::with_progress({
    p <- progressr::progressor(steps = length(row_combs))
    future::plan('future::multisession', workers = n)
    res <- row_combs |> furrr::future_map(\(x) {
        p()
        i <- x[1]; j <- x[2]
        spa <- betaPDF_ps(par_space[i, ])
        spb <- betaPDF_ps(par_space[j, ])
        ov <- beta_overlap(spa, spb)
        results[1, ] <- c(i, j,
                          niche_diss(spa, spb, ov),
                          niche_excl(spa, spb),
                          par_space[i, ], par_space[j, ])
        results
    }) |>
        dplyr::bind_rows()
    future::plan(future::sequential)
    })
    res
    }


#' Niche Divergence Calculations for a given Parameter Space based on the beta function curve
#' @description
#' This functions calculates niche divergence indices on all theoretical response curves simulated using the set_parameter_space(). It is a wrapper around two other functions that work in single or parallel processing.
#'
#' @param psp Parameter space created using the set_parameter_space() function
#' @param parallel logical. Change to TRUE if parallel processing is desired
#' @param n number of cores to use in parallel processing (see future::plan())
#'
#' @return A data frame containing niche divergence plane values of Dissimilarity and Exclusivity, along with the parameters of each species response
#' @export
#'
#' @examples
#' # Set parameter domains
#' a = b <- round(1:3, 3)
#' alpha = gamma <- round(1)
#'
#' # Set parameter space (combination of parameters for any beta function)
#' par_space <- set_parameter_space(a, b, alpha, gamma)
#'
#' #Run calc
#' calculate_nichediv_beta(psp = par_space) # single core calculation
#' calculate_nichediv_beta(psp = par_space, parallel = TRUE, n = 2) # parallel calculations
#' #By default function uses 2 cores, but you can run
#' #n = parallelly::availableCores()-1
calculate_nichediv_beta <- function(psp, parallel = FALSE, n = 2) {
    if(parallel) {
        calculate_nichediv_bpar(psp = psp, n = n)
    } else {
        calculate_nichediv_b1(psp = psp)
    }
}
