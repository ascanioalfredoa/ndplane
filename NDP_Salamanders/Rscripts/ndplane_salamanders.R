#' Load Salamander Maps Data
#'
#' Loads spatial vector data for salamander occurrences.
#' @param filepath Character. File path to the shapefile.
#' @return Spatial vector object.
#' @import terra
#' @export
load_salamander_maps <- function(filepath) {
  terra::vect(filepath)
}

#' Load Salamander CSV Data
#'
#' Loads salamander environmental data CSV.
#' @param filepath Character. File path to the CSV file.
#' @return Data frame.
#' @import readr
#' @export
load_salamander_data <- function(filepath) {
  readr::read_csv(filepath)
}

#' Calculate Density Curves
#'
#' Calculate density curves for each species category over selected variables.
#' @param dat Data frame with environmental data.
#' @param category_col Character. Name of the column for species categories.
#' @param variables Character vector. Variables to analyze.
#' @return List with densities for each variable and category.
#' @import dplyr
#' @import purrr
#' @export
calculate_density_curves <- function(dat, category_col, variables) {
  densities <- list()
  for (var in variables) {
    envdat <- dat[!is.na(dat[[var]]), ]
    x_all <- density(envdat[[var]], n = 3000)$x
    cat_names <- envdat |> dplyr::select({{category_col}}) |> unique() |> dplyr::pull()
    ndp_data <- envdat |>
      dplyr::select({{category_col}}, {{var}}) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(category_col))) |>
      dplyr::group_map(~ {
        dens <- density(.x[[var]], n = 3000)
        data.frame(x = dens$x, y = dens$y)
      }) |>
      rlang::set_names(cat_names) |>
      lapply(function(x) {
        stats::approx(x$x, x$y, xout = x_all) |>
          as.data.frame() |>
          dplyr::filter(complete.cases(.)) |>
          dplyr::mutate(Cumden = cumsum(y)) |>
          dplyr::mutate(Cumden = Cumden / max(Cumden)) |>
          dplyr::filter(Cumden >= 0.025 & Cumden <= 0.975) |>
          dplyr::mutate(y = y / max(y)) |>
          dplyr::select(1:2)
      }) |>
      dplyr::bind_rows(.id = category_col)
    densities[[var]] <- ndp_data
  }
  densities
}

#' Calculate Niche Dissimilarity Metrics
#'
#' Calculate niche dissimilarity (NicheDiss) and niche exclusivity (NicheEx) between category pairs per variable.
#' @param densities List of density data frames per variable with category column.
#' @param category_col Character. Name of category column.
#' @return Data frame with NicheDiss and NicheEx.
#' @import dplyr
#' @import purrr
#' @importFrom pracma trapz
#' @export
calculate_niche_metrics <- function(densities, category_col) {
  res <- tibble::tibble(Taxa1 = character(), Taxa2 = character(), Var = character(), NicheDiss = double(), NicheEx = double())
  for (var in names(densities)) {
    ndp_data <- densities[[var]]
    taxa_comb <- t(combn(unique(ndp_data[[category_col]]), 2))
    for (i in seq_len(nrow(taxa_comb))) {
      A <- ndp_data |> dplyr::filter(!!rlang::sym(category_col) == taxa_comb[i, 1])
      B <- ndp_data |> dplyr::filter(!!rlang::sym(category_col) == taxa_comb[i, 2])
      C <- dplyr::inner_join(A, B, by = "x", suffix = c("_a", "_b")) |>
        dplyr::transmute(x, y = pmin(y_a, y_b))
      NicheDiss <- 1 - (pracma::trapz(C$x, C$y) / pracma::trapz(A$x, A$y) +
        pracma::trapz(C$x, C$y) / pracma::trapz(B$x, B$y)) / 2
      RNO <- 1 - (min(max(A$x), max(B$x)) - max(min(A$x), min(B$x))) /
        (max(max(A$x), max(B$x)) - min(min(A$x), min(B$x)))
      if (RNO > 1) RNO <- 1
      res <- dplyr::bind_rows(res, tibble::tibble(
        Taxa1 = taxa_comb[i, 1],
        Taxa2 = taxa_comb[i, 2],
        Var = var,
        NicheDiss = NicheDiss,
        NicheEx = RNO
      ))
    }
  }
  res
}

#' Plot Niche Divergence Plane
#'
#' Create a ggplot object visualizing niche exclusivity vs dissimilarity.
#' @param niche_metrics Data frame from calculate_niche_metrics().
#' @return ggplot2 plot object.
#' @import ggplot2
#' @export
plot_ndp <- function(niche_metrics) {
  niche_metrics <- niche_metrics |>
    dplyr::mutate(comparison = paste(Taxa1, Taxa2, sep = "/")) |>
    dplyr::mutate(vartype = c(rep("Temperature", 11), rep("Precipitation", 8), rep("Soil", 12)))
  ggplot2::ggplot(niche_metrics, aes(x = NicheEx, y = NicheDiss, label = Var, fill = vartype)) +
    ggplot2::xlab("Niche Exclusivity") + ggplot2::ylab("Niche Dissimilarity") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::geom_hline(yintercept = 0.5) +
    ggplot2::geom_vline(xintercept = 0.5) +
    ggplot2::geom_point(size = 4, shape = 21, alpha = 0.8) +
    ggrepel::geom_text_repel(size = 4) +
    ggplot2::scale_fill_manual("Variable type", values = c("#62A39F", "#000000", "#2F8745")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = c(0.87, 0.75), legend.background = ggplot2::element_rect(fill = "white", color = "black"))
}

#' Save Plot to Files
#'
#' Save a given plot to PDF and PNG.
#' @param plot_obj ggplot object.
#' @param pdf_file Path for PDF output.
#' @param png_file Path for PNG output.
#' @param width Width for PNG in pixels.
#' @param height Height for PNG in pixels.
#' @param res Resolution for PNG.
#' @export
save_plots <- function(plot_obj, pdf_file, png_file, width = 1900, height = 1100, res = 300) {
  ggplot2::ggsave(pdf_file, plot = plot_obj, device = "pdf", width = width / 100, height = height / 100, units = "in")
  png(png_file, width = width, height = height, res = res)
  print(plot_obj)
  dev.off()
}

#' Full Analysis Pipeline for Salamander Niche Divergence Plane
#'
#' Runs full pipeline from loading, calculating densities, niche metrics and plotting.
#' @param shapefile_path Path to salamander occurrences shapefile.
#' @param csv_data_path Path to salamander environmental CSV data.
#' @param output_pdf Path to save PDF plot.
#' @param output_png Path to save PNG plot.
#' @param category_col Name of species category column.
#' @param variables Variables to analyze. If NULL, auto-detects.
#' @return List with densities, niche metrics, and plot object.
#' @export
run_sal_ndp_analysis <- function(shapefile_path, csv_data_path, output_pdf, output_png, category_col = "Cmmn_Nm", variables = NULL) {
  maps <- load_salamander_maps(shapefile_path)
  dat <- load_salamander_data(csv_data_path)
  if (is.null(variables)) {
    variables <- names(dat)[-(1:2)]
    variables <- variables[!grepl("cl$|_Levels", variables)]
  }
  densities <- calculate_density_curves(dat, category_col, variables)
  niche_metrics <- calculate_niche_metrics(densities, category_col)
  plot_obj <- plot_ndp(niche_metrics)
  save_plots(plot_obj, output_pdf, output_png)
  list(densities = densities, niche_metrics = niche_metrics, plot = plot_obj)
}
