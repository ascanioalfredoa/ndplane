library(tictoc)

# Set parameter domains
a = b <- round(1:3, 3)
alpha = gamma <- round(1:3)
#alpha = gamma <- seq(0.1, 10, 0.5)

# Set parameter space (combination of parameters for any beta function)
par_space <- set_parameter_space(a, b, alpha, gamma)

#### Apply beta functions to parameter space ####
row_combs <- combn(1:nrow(par_space), 2, simplify = F)

results_names <- c("spA_par", "spB_par", "Dissimilarity", "Exclusivity",
                   paste(names(par_space), "_A", sep = ""),
                   paste(names(par_space), "_B", sep = "")
)

results <- dplyr::as_tibble(array(data = 0,
                               dim = c(1, length(results_names)),
                               dimnames = list(1, results_names)))


tic()
progressr::with_progress({
    p <- progressr::progressor(steps = length(row_combs))
    row_combs |> purrr::map_df(\(x) {
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
    })
}
)
toc()

tic()
p <- progressor(steps = nrow(row_combs))
row_combs |> purrr::map_df(\(x) {
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
})
toc()


tic()
progressr::with_progress({
    p <- progressr::progressor(steps = length(row_combs))
future::plan(future::multisession, workers = 10)
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
}, .progress = TRUE) |>
    dplyr::bind_rows()
future::plan(future::sequential)
})
toc()
