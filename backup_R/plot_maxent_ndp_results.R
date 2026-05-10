#' Plot Niche Divergence Plane indices
#'
#' @param ind_res Data frame containing columns dissimilarity, exclusivity, and var
#' @param exclude_vars Character vector of variable names to exclude from the plot (default c("bio8", "bio9"))
#' @return ggplot2 object with niche exclusivity vs dissimilarity scatterplot
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline scale_fill_manual scale_shape_manual theme_bw theme labs xlim ylim
#' @importFrom ggrepel geom_text_repel
#' @export
plot_ndp_indices <- function(ind_res, exclude_vars = c("bio8", "bio9")) {

  classify_var_type <- function(vars) {
    ifelse(vars %in% paste0("bio", 1:11), "Temperature",
           ifelse(vars %in% paste0("bio", 12:19), "Precipitation", "Soil"))
  }

  ind_res <- ind_res |>
    dplyr::filter(!var %in% exclude_vars) |>
    dplyr::mutate(vartype = classify_var_type(var))

  p <- ggplot2::ggplot(ind_res, aes(x = exclusivity, y = dissimilarity, label = var, fill = vartype, shape = vartype)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed") +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggrepel::geom_text_repel(size = 3) +
    ggplot2::scale_fill_manual("Variable type", values = c("Temperature" = "#62A39F",
                                                 "Precipitation" = "#000000",
                                                 "Soil" = "#2F8745")) +
    ggplot2::scale_shape_manual("Variable type", values = c(21, 22, 23)) +
    ggplot2::labs(x = "Niche Exclusivity", y = "Niche Dissimilarity") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = c(0.87, 0.75),
          legend.background = element_rect(fill = "white", color = "black"))

  return(p)
}

#' Plot Maxent response curves
#'
#' @param plot_res Data frame containing response curve data with columns layer, value, source, var
#' @param selected_vars Character vector of variables to plot (default c("bio1", "bio3", "bio12", "bio14", "Sand", "WC3rdbar"))
#' @param color_palette Named vector of species colors (default c("#000000", "#62A39F"))
#' @return ggplot2 object with faceted line plots of response curves by variable and type
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap ylab xlab scale_color_manual theme_bw theme
#' @importFrom dplyr mutate filter .data
#' @importFrom forcats fct_recode fct_relevel
#' @export
plot_response_curves <- function(plot_res,
                                 selected_vars = c("bio1", "bio3", "bio12", "bio14", "Sand", "WC3rdbar"),
                                 color_palette = c("#000000", "#62A39F")) {

  classify_var_type <- function(vars) {
    ifelse(vars %in% paste0("bio", 1:11), "Temperature",
           ifelse(vars %in% paste0("bio", 12:19), "Precipitation", "Soil"))
  }

  var_labels <- c(
    bio1 = "Annual Mean Temperature (\u00B0 C) [bio1]",
    bio3 = "Isothermality (\u00B0 C) [bio3]",
    bio12 = "Annual Precipitation (kg/m^2) [bio12]",
    bio14 = "Precipitation of Driest Month (kg/m^2) [bio14]",
    Sand = "Sand Total (%)",
    WC3rdbar = "Water Content at 1/3 Bar (%) [WC3rdbar]"
  )

  plot_res2 <- plot_res |>
    dplyr::filter(.data$var %in% selected_vars) |>
    dplyr::mutate(
      var = factor(.data$var, levels = selected_vars),
      var = forcats::fct_relevel(.data$var, selected_vars),
      var = forcats::fct_recode(.data$var, !!!var_labels),
      vartype = classify_var_type(as.character(.data$var))
    )

  p <- ggplot2::ggplot(plot_res2, aes(x = layer, y = value, color = source)) +
    ggplot2::geom_line(linewidth = 1.2, alpha = 0.8) +
    ggplot2::facet_wrap(vartype ~ var, scales = "free_x") +
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("Environmental gradient") +
    ggplot2::scale_color_manual("Species", values = color_palette) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")

  return(p)
}

#' Wrap Niche Divergence Plane plotting: indices and response curves, optionally combined with species occurrence map
#'
#' @param ind_res Data frame with niche indices per variable
#' @param plot_res Data frame of combined response curves
#' @param occ_sf Optional sf/spatial object with species occurrences for map plotting
#' @param species_col Optional character, name of species column in occ_sf for map color
#' @param base_map_sf Optional sf polygon layer for basemap in species map plot
#' @param selected_vars Character vector of variables to highlight in response curves plot
#' @param color_palette Named vector of species colors
#' @return named list with ggplot objects: map_plot, ndp_plot, response_curves_plot, combined_plot (if map provided)
#' @importFrom cowplot draw_plot_label
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggpubr as_ggplot
#' @export
plot_maxent_ndp_all <- function(ind_res, plot_res,
                               occ_sf = NULL, species_col = NULL, base_map_sf = NULL,
                               selected_vars = c("bio1", "bio3", "bio12", "bio14", "Sand", "WC3rdbar"),
                               color_palette = c("#000000", "#62A39F")) {

  p_ndp <- plot_ndp_indices(ind_res)
  p_resp <- plot_response_curves(plot_res, selected_vars, color_palette)

  p_map <- NULL
  if (!is.null(occ_sf) && !is.null(base_map_sf)) {
    if (!species_col %in% names(occ_sf)) {
      warning("species_col not found in occ_sf; skipping map plot")
    } else {
      p_map <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = base_map_sf, fill = "grey90", color = NA) +
        ggplot2::geom_sf(data = occ_sf, ggplot2::aes_string(color = species_col), alpha = 0.5, size = 0.5) +
        ggplot2::geom_sf(data = base_map_sf, fill = NA, color = "grey20", size = 0.3) +
        ggplot2::scale_color_manual(values = color_palette, name = "Species") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1)))
    }
  }

  combined_plot <- NULL
  if (!is.null(p_map)) {
    combined <- gridExtra::arrangeGrob(
      p_map + ggplot2::theme(axis.title = ggplot2::element_text(size = 8),
                    axis.text = ggplot2::element_text(size = 6),
                    legend.text = ggplot2::element_text(size = 6),
                    legend.title = ggplot2::element_text(size = 8),
                    legend.spacing.x = ggplot2::unit(0.15, "cm"),
                    legend.key.size = ggplot2::unit(0.2, "cm")) +
        ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 1))),
      p_ndp + ggplot2::theme(axis.title = ggplot2::element_text(size = 8),
                    legend.text = ggplot2::element_text(size = 8),
                    text = ggplot2::element_text(size = 8)),
      p_resp + ggplot2::theme(text = ggplot2::element_text(size = 8),
                     legend.position = "none",
                     axis.title = ggplot2::element_text(size = 8)),
      ncol = 2, nrow = 2,
      layout_matrix = rbind(c(1,1,3,3,3,3,3),
                            c(1,1,3,3,3,3,3),
                            c(1,1,3,3,3,3,3),
                            c(2,2,2,2,2,2,2),
                            c(2,2,2,2,2,2,2),
                            c(2,2,2,2,2,2,2),
                            c(2,2,2,2,2,2,2),
                            c(2,2,2,2,2,2,2))
    )
    combined_plot <- ggpubr::as_ggplot(combined) +
      cowplot::draw_plot_label(label = c("A)", "B)", "C)"), size = 10,
                              x = c(-0.005, 0.30, -0.005), y = c(1, 1, 0.62))
  }

  return(list(
    map_plot = p_map,
    ndp_plot = p_ndp,
    response_curves_plot = p_resp,
    combined_plot = combined_plot
  ))
}