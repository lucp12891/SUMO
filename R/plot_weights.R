#' @name plot_weights
#' @title Visualize feature loadings (weights)
#' @description Generate scatter or histogram plots of feature loadings (weights)
#'   from simulated or real multi-omics data. Supports per-omic views and, when
#'   available, an integrated view.
#'
#' @param sim_object A multi-omics object (e.g., from \code{simulate_twoOmicsData()}
#'   and \code{as_multiomics()}).
#' @param omic Integer or character. Which view to plot: \code{1} (omic.one),
#'   \code{2} (omic.two), or \code{"integrated"} (if present). Default \code{1}.
#' @param factor_num Integer or "all". Which factor(s) to visualize. Default \code{1}.
#' @param type Character. Plot type: \code{"scatter"} or \code{"histogram"}.
#'   Default \code{"scatter"}.
#' @param show.legend Logical. Whether to show the legend. Default \code{TRUE}.
#'
#' @return A \link[ggplot2]{ggplot} object (single plot) or a grob returned by
#'   \link[gridExtra]{grid.arrange} when multiple panels are combined.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point after_stat labs
#' @importFrom ggplot2 scale_color_viridis_c scale_fill_viridis_c theme
#' @importFrom ggplot2 element_text element_blank theme_bw theme_minimal ggtitle
#' @importFrom gridExtra grid.arrange
#' @importFrom rlang .data
#'
#' @examples
#' output_obj <- simulate_twoOmicsData(
#'   vector_features = c(4000, 3000),
#'   n_samples = 100,
#'   n_factors = 2,
#'   signal.samples = NULL,
#'   signal.features.one = NULL,
#'   signal.features.two = NULL,
#'   snr = 2.5,
#'   num.factor = "multiple",
#'   advanced_dist = "mixed"
#' )
#'
#' output_obj <- as_multiomics(output_obj)
#'
#' plot_weights(
#'   sim_object = output_obj,
#'   factor_num = 1,
#'   omic = 1,
#'   type = "scatter",
#'   show.legend = FALSE
#' )
#'
#' plot_weights(
#'   sim_object = output_obj,
#'   factor_num = 2,
#'   omic = 2,
#'   type = "histogram"
#' )
#'
#' @export
plot_weights <- function(sim_object,
                         omic = 1,
                         factor_num = 1,
                         type = "scatter",
                         show.legend = TRUE) {

  # --- helpers to detect structure ------------------------------------------------
  .is_multi_omics <- function(x) {
    is.list(x$list_betas) && length(x$list_betas) >= 1 &&
      all(vapply(x$list_betas, is.list, logical(1)))
  }

  # Normalize omic index and build feature name prefix
  if (is.character(omic)) {
    omic_idx <- as.integer(gsub("[^0-9]", "", omic))
  } else {
    omic_idx <- as.integer(omic)
  }
  if (is.na(omic_idx) || omic_idx < 1) stop("`omic` must be a positive integer (e.g., 1, 2) or like 'omic1'.")

  feature_prefix <- paste0("^omic", omic_idx, "_feature_")

  # --- pull features from concatenated data (works for multi-omics) --------------
  if (!is.null(sim_object$concatenated_datasets) && length(sim_object$concatenated_datasets) >= 1) {
    sim_data <- as.data.frame(sim_object$concatenated_datasets[[1]])
    features <- grep(feature_prefix, colnames(sim_data), value = TRUE)
  } else if (!is.null(sim_object$sim_data_example)) {
    # very old structure fallback
    sim_data <- as.data.frame(sim_object$sim_data_example)
    features <- grep(feature_prefix, colnames(sim_data), value = TRUE)
  } else {
    stop("Could not find concatenated feature matrix in `sim_object`.")
  }

  # --- pick the right betas list -------------------------------------------------
  if (.is_multi_omics(sim_object)) {
    # simulateMultiOmics(): list_betas is a list per omic, each containing beta1, beta2, ...
    if (omic_idx > length(sim_object$list_betas))
      stop("Requested `omic` not available in `list_betas`.")
    factor_list <- sim_object$list_betas[[omic_idx]]
  } else {
    # Backward compatibility: two-omics object with list_betas (omic1) and list_deltas (omic2)
    if (omic_idx == 1) {
      factor_list <- sim_object$list_betas
    } else if (omic_idx == 2 && !is.null(sim_object$list_deltas)) {
      factor_list <- sim_object$list_deltas
    } else {
      stop("For two-omics objects, use omic = 1 (betas) or omic = 2 (deltas).")
    }
  }

  if (is.null(factor_list) || length(factor_list) == 0)
    stop("No factor loadings found for the selected omic.")

  # factor ids from names like 'beta1', 'beta2', ...
  fl_names <- names(factor_list)
  if (is.null(fl_names)) fl_names <- rep("", length(factor_list))
  factor_ids <- suppressWarnings(as.integer(gsub("[^0-9]", "", fl_names)))
  # if names missing, assume sequential 1..k
  if (any(is.na(factor_ids))) factor_ids <- seq_along(factor_list)
  names(factor_list) <- paste0("beta", factor_ids)

  # ---- plotting core ------------------------------------------------------------
  build_df <- function(fid) {
    vec <- factor_list[[paste0("beta", fid)]]
    if (is.null(vec)) stop(paste0("Factor ", fid, " not found for selected omic."))
    df <- data.frame(features = features[seq_along(vec)], factor = as.numeric(vec))
    df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", df$features))
    df[order(df$Index), ]
  }

  make_one_plot <- function(fid) {
    df <- build_df(fid)
    ttl <- paste0("Omic ", omic_idx, " | Factor ", fid)
    if (identical(type, "scatter")) {
      build_scatter_plot(df, ttl, show.legend)
    } else if (identical(type, "histogram")) {
      build_histogram_plot(df, ttl, show.legend)
    } else {
      stop("`type` must be 'scatter' or 'histogram'.")
    }
  }

  if (identical(factor_num, "all")) {
    plot_list <- lapply(factor_ids, make_one_plot)
    names(plot_list) <- paste0("F", factor_ids)
    return(do.call(gridExtra::grid.arrange, c(plot_list, ncol = min(3, length(plot_list)))))
  }

  if (!is.numeric(factor_num) || length(factor_num) != 1)
    stop("`factor_num` must be a single integer or 'all'.")
  if (!(factor_num %in% factor_ids))
    stop("Requested factor not available for this omic.")

  make_one_plot(as.integer(factor_num))
}

#' Helper function to build scatter plot (internal use)
#' @param df Dataframe with features and weights
#' @param title Plot title
#' @param show.legend Logical to show legend
#' @keywords internal
build_scatter_plot <- function(df, title = "", show.legend = TRUE) {
  p <- ggplot(df, aes(x = Index, y = factor, color = factor)) +
    geom_point() +
    labs(x = "Features", y = "Weights", color = "") +
    scale_color_viridis_c() +
    theme_bw() +
    ggtitle(title) +
    theme(plot.title = element_text(size = 10, hjust = 0.5))

  if (!show.legend) p <- p + theme(legend.position = "none")
  return(p)
}

#' Helper function to build histogram plot (internal use)
#' @param df Dataframe with features and weights
#' @param title Plot title
#' @param show.legend Logical to show legend
#' @keywords internal
build_histogram_plot <- function(df, title = "", show.legend = TRUE) {
  p <- ggplot(df, aes(x = factor, fill = after_stat(x))) +
    geom_histogram(bins = 100) +
    labs(x = "Loadings", y = "Frequency", fill = "") +
    scale_fill_viridis_c() +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = if (show.legend) "right" else "none",
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle(title)
  return(p)
}

utils::globalVariables(c("Index", "x"))
# # Factor 1 of omic 2, scatter
# plot_weights(sim_data_example, omic = 2, factor_num = 2, type = "scatter", show.legend = FALSE)
#
# # All factors of omic 1, histograms
# plot_weights(sim_data_example, omic = 2, factor_num = 1, type = "histogram")
