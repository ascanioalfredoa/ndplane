test_that("normalize_suitability returns normalized raster", {
    skip_if_not_installed("virtualspecies")
    skip_if_not_installed("terra")

    # Simulate basic environment and species
    r <- terra::rast(nrows = 10, ncols = 10, vals = runif(100))
    names(r) <- "bio5"

    params = virtualspecies::formatFunctions(bio5 = c(fun = "dnorm", mean = 0.5, sd = 0.2))

    vsp <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = params)

    norm <- normalize_suitability(vsp)
    expect_true(inherits(norm, "SpatRaster"))
    expect_true(abs(sum(terra::values(norm), na.rm = TRUE) - 1) < 0.01)
})

test_that("calculate_similarity computes D and I correctly", {
    skip_if_not_installed("terra")

    r <- terra::rast(nrows = 10, ncols = 10, vals = runif(100))
    names(r) <- "bio5"

    params1 = virtualspecies::formatFunctions(bio5 = c(fun = "dnorm", mean = 0.4, sd = 0.2))
    params2 = virtualspecies::formatFunctions(bio5 = c(fun = "dnorm", mean = 0.6, sd = 0.2))

    vsp1 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = params1) |> normalize_suitability()
    vsp2 <- virtualspecies::generateSpFromFun(raster.stack = r, parameters = params2) |> normalize_suitability()

    sim <- calculate_similarity(vsp1, vsp2, method = c("D", "I"))
    expect_named(sim, c("D", "I"))
    expect_true(all(sim <= 1))
    expect_true(all(sim >= 0))
})

test_that("define_vsp_param_space creates valid rasters and parameter grids", {
    param_ranges <- list(
        temp = list(a = seq(10, 15, 5), b = seq(20, 25, 5), alpha = c(1, 2), gamma = c(1, 2)),
        precip = list(a = seq(500, 600, 100), b = seq(700, 800, 100), alpha = c(2, 3), gamma = c(2, 3))
    )

    res <- define_vsp_param_space(param_ranges)

    # Structure
    expect_named(res, c("env_raster", "parameter_grids"))
    expect_type(res$parameter_grids, "list")
    expect_s4_class(res$env_raster, "SpatRaster")

    # Check raster layer names
    expect_equal(sort(names(res$env_raster)), sort(c("temp", "precip")))
    expect_equal(terra::nlyr(res$env_raster), 2)

    # Check parameter grids
    for (grid in res$parameter_grids) {
        expect_s3_class(grid, "tbl_df")
        expect_true(all(grid$a < grid$b))
    }
})

test_that("define_vsp_param_space accepts a custom raster stack", {
    r1 <- terra::rast(nrows = 10, ncols = 10, vals = runif(100))
    r2 <- terra::rast(nrows = 10, ncols = 10, vals = runif(100))
    names(r1) <- "temp"
    names(r2) <- "precip"
    real_stack <- c(r1, r2)

    param_ranges <- list(
        temp = list(a = seq(10, 15, 5), b = seq(20, 25, 5), alpha = c(1), gamma = c(1)),
        precip = list(a = seq(500, 600, 100), b = seq(700, 800, 100), alpha = c(2), gamma = c(2))
    )

    res <- define_vsp_param_space(param_ranges, raster_stack = real_stack)
    expect_equal(terra::nlyr(res$env_raster), 2)
    expect_equal(dim(res$env_raster)[1:2], c(10, 10))
})

test_that("define_vsp_param_space writes parameter grids to files", {
    tempdir_out <- tempfile()
    dir.create(tempdir_out)

    param_ranges <- list(
        temp = list(a = seq(10, 15, 5), b = seq(20, 25, 5), alpha = c(1), gamma = c(1)),
        precip = list(a = seq(500, 600, 100), b = seq(700, 800, 100), alpha = c(2), gamma = c(2))
    )

    res <- define_vsp_param_space(param_ranges, save_to = tempdir_out)

    temp_file <- file.path(tempdir_out, "params_temp.csv")
    precip_file <- file.path(tempdir_out, "params_precip.csv")

    expect_true(file.exists(temp_file))
    expect_true(file.exists(precip_file))

    df <- readr::read_csv(temp_file, show_col_types = FALSE)
    expect_named(df, c("a", "b", "alpha", "gamma"))
    unlink(tempdir_out, recursive = TRUE)
})

test_that("sample_params returns one row per variable", {
    param_grids <- list(
        temp = tibble::tibble(a = 1:5, b = 6:10, alpha = 1, gamma = 1),
        precip = tibble::tibble(a = 10:14, b = 15:19, alpha = 2, gamma = 2)
    )

    sampled <- sample_params(param_grids)
    expect_named(sampled, c("temp", "precip"))
    expect_true(all(purrr::map_int(sampled, nrow) == 1))
})

test_that("format_vsp returns a valid virtualspecies list", {
    sampled <- list(
        temp = tibble::tibble(a = 10, b = 20, alpha = 2, gamma = 2),
        precip = tibble::tibble(a = 100, b = 300, alpha = 3, gamma = 3)
    )

    out <- format_vsp(sampled)
    expect_type(out, "list")
    expect_named(out, c("temp", "precip"))
    expect_true(all(sapply(out, function(x) "fun" %in% names(x))))
})

test_that("nd_indices returns all NDP values", {
    pA <- tibble::tibble(a = 10, b = 20, alpha = 2, gamma = 2)
    pB <- tibble::tibble(a = 15, b = 25, alpha = 3, gamma = 2)

    out <- nd_indices(pA, pB)
    expect_named(out, c("ND", "NE", "ND_magnitude", "ND_angle"))
    expect_s3_class(out, "tbl_df")
    expect_true(all(out$ND >= 0 & out$NE >= 0))
})

test_that("simulate_vsp generates species and returns D/I", {
    r1 <- terra::rast(nrows = 10, ncols = 10, vals = runif(100, 10, 30))
    r2 <- terra::rast(nrows = 10, ncols = 10, vals = runif(100, 500, 2000))
    names(r1) <- "temp"
    names(r2) <- "precip"
    env <- c(r1, r2)

    params1 <- format_vsp(list(
        temp = tibble::tibble(a = 10, b = 30, alpha = 2, gamma = 2),
        precip = tibble::tibble(a = 500, b = 2000, alpha = 2, gamma = 2)
    ))
    params2 <- format_vsp(list(
        temp = tibble::tibble(a = 15, b = 35, alpha = 3, gamma = 3),
        precip = tibble::tibble(a = 600, b = 1800, alpha = 3, gamma = 3)
    ))

    sim <- simulate_vsp(params1, params2, env)
    expect_true(sim$D >= 0 && sim$D <= 1)
    expect_true(sim$I >= 0 && sim$I <= 1)
    expect_s4_class(sim$sp1, "SpatRaster")
    expect_s4_class(sim$sp2, "SpatRaster")
})

test_that("run_vsp_sim produces valid simulation output", {
    param_ranges <- list(
        temp = list(a = seq(10, 20, 5), b = seq(25, 30, 5), alpha = c(2), gamma = c(2)),
        precip = list(a = seq(500, 800, 100), b = seq(1000, 2000, 500), alpha = c(2), gamma = c(2))
    )

    vsp <- define_vsp_param_space(param_ranges)
    sim_results <- run_vsp_sim(vsp, n_comparisons = 3)

    expect_s3_class(sim_results, "tbl_df")
    expect_true(all(c("comparison_id", "variable", "ND", "NE", "D", "I") %in% names(sim_results)))
    expect_equal(length(unique(sim_results$comparison_id)), 3)
    expect_true(all(sim_results$ND >= 0 & sim_results$ND <= 1))
    expect_true(all(sim_results$I >= 0 & sim_results$I <= 1))
})

test_that("run_vsp_sim respects seed and returns reproducible results", {
    param_ranges <- list(
        temp = list(a = seq(10, 20, 5), b = seq(25, 30, 5), alpha = c(2), gamma = c(2)),
        precip = list(a = seq(500, 800, 100), b = seq(1000, 2000, 500), alpha = c(2), gamma = c(2))
    )

    vsp <- define_vsp_param_space(param_ranges)
    sim1 <- run_vsp_sim(vsp, n_comparisons = 3, seed = 42, progress = FALSE)
    sim2 <- run_vsp_sim(vsp, n_comparisons = 3, seed = 42, progress = FALSE)

    expect_equal(sim1, sim2)
})

test_that("run_vsp_sim can export CSV and RDS files", {
    param_ranges <- list(
        temp = list(a = seq(10, 20, 5), b = seq(25, 30, 5), alpha = c(2), gamma = c(2)),
        precip = list(a = seq(500, 800, 100), b = seq(1000, 2000, 500), alpha = c(2), gamma = c(2))
    )

    vsp <- define_vsp_param_space(param_ranges)

    tmp_csv <- tempfile(fileext = ".csv")
    tmp_rds <- tempfile(fileext = ".rds")

    run_vsp_sim(vsp, n_comparisons = 2, save_as = tmp_csv, progress = FALSE)
    run_vsp_sim(vsp, n_comparisons = 2, save_as = tmp_rds, progress = FALSE)

    expect_true(file.exists(tmp_csv))
    expect_true(file.exists(tmp_rds))

    df_csv <- readr::read_csv(tmp_csv, show_col_types = FALSE)
    df_rds <- readRDS(tmp_rds)

    expect_s3_class(df_csv, "tbl_df")
    expect_s3_class(df_rds, "tbl_df")

    unlink(c(tmp_csv, tmp_rds))
})
