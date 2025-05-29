#' @name plot_weights
#' @title Visualizing the raw loading/weights of the features
#' @description Generates scatter or histogram plots of feature loadings (weights) from simulated multi-omics data.
#' Supports plotting for omic.one, omic.two, or integrated views.
#' @param sim_object R object containing data to be plotted.
#' @param factor_num Integer or "all". Specifies which factor(s) to visualize.
#' @param data Character. Section of the data to visualize: "omic.one", "omic.two", or "integrated".
#' @param type Character. Type of plot: "scatter" or "histogram".
#' @param show.legend Logical. Whether to show the legend in the plot. Default is TRUE.
#' @return A ggplot object or a combined grid of plots.
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point after_stat labs scale_color_viridis_c scale_fill_viridis_c theme element_text element_blank theme_bw theme_minimal ggtitle element_text element_blank
#' @importFrom gridExtra grid.arrange
#' @importFrom readxl read_excel
#' @importFrom readr read_csv
#' @importFrom stringr str_detect
#' @importFrom rlang .data
#' @examples
#' output_obj <- simulate_twoOmicsData(
#'   vector_features = c(4000, 3000),
#'   n_samples = 100,
#'   n_factors = 2,
#'   signal.samples = NULL,
#'   signal.features.one = NULL,
#'   signal.features.two = NULL,
#'   snr = 2.5,
#'   num.factor = 'multiple',
#'   advanced_dist = 'mixed'
#'   )
#'
#' plot_weights(
#'    sim_object = output_obj,
#'    factor_num = 1,
#'    data = 'omic.one',
#'    type = 'scatter',
#'    show.legend = FALSE
#'    )
#'
#' plot_weights(
#'   sim_object = output_obj,
#'   factor_num = 2,
#'   data = 'omic.two',
#'   type = 'histogram'
#'   )
#' @export
plot_weights <- function(sim_object, factor_num = 1, data = 'omic.one', type = 'scatter', show.legend = TRUE) {
  sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

  if (data == 'omic.one') {
    features <- grep("^omic1_feature_", colnames(sim_data), value = TRUE)
    factor_list <- sim_object$list_betas
  } else if (data == 'omic.two') {
    features <- grep("^omic2_feature_", colnames(sim_data), value = TRUE)
    factor_list <- sim_object$list_deltas
  } else {
    stop("Invalid data section. Choose from 'omic.one' or 'omic.two'.")
  }

  factor_ids <- as.numeric(gsub("[^0-9]", "", names(factor_list)))
  plot_list <- list()

  if (factor_num == "all") {
    for (i in factor_ids) {
      factor <- factor_list[[paste0(ifelse(data == 'omic.one', 'beta', 'delta'), i)]]
      df <- data.frame(features = features, factor = factor)
      df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", df$features))
      df <- df[order(df$Index), ]
      title <- paste("Data:", data, "| Factor:", i)

      p <- if (type == 'scatter') {
        build_scatter_plot(df, title, show.legend)
      } else {
        build_histogram_plot(df, title, show.legend)
      }
      plot_list[[i]] <- p
    }
    #return(gridExtra::grid.arrange(grobs = plot_list, ncol = length(plot_list)))
    return(do.call(gridExtra::grid.arrange, c(plot_list, ncol = length(plot_list))))

  } else if (factor_num %in% factor_ids) {
    i <- factor_num
    factor <- factor_list[[paste0(ifelse(data == 'omic.one', 'beta', 'delta'), i)]]
    df <- data.frame(features = features, factor = factor)
    df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", df$features))
    df <- df[order(df$Index), ]
    title <- paste("Data:", data, "| Factor:", i)

    if (type == 'scatter') {
      return(build_scatter_plot(df, title, show.legend))
    } else {
      return(build_histogram_plot(df, title, show.legend))
    }
  } else {
    stop("Invalid factor_num. Provide an integer matching available factor numbers or 'all'.")
  }
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
