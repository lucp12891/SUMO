#' @name plot_factor
#' @title Visualization of factor scores  (ground truth)
#' @description Scatter or histogram plots of sample-level factor scores from simulated multi-omics data, using scores from list_alphas and list_gammas.
#' @param sim_object R object containing simulated data output from OmixCraftHD
#' @param factor_num Integer or "all". Which factor(s) to plot.
#' @param type Character. Either "scatter" (default) or "histogram" for plot type.
#' @param show.legend Logical. Whether to show legend in plots. Default is TRUE.
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point after_stat labs scale_color_viridis_c scale_fill_viridis_c theme element_text element_blank theme_bw theme_minimal ggtitle element_text element_blank
#' @importFrom gridExtra grid.arrange
#' @importFrom readxl read_excel
#' @importFrom readr read_csv
#' @importFrom stringr str_detect
#' @examples
#' output_obj <- OmixCraftHD(
#'   vector_features = c(4000,3000),
#'   n_samples = 100,
#'   n_factors = 2,
#'   snr = 2.5,
#'   num.factor = 'multiple',
#'   advanced_dist = 'mixed')
#'
#' plot_factor(sim_object = output_obj, factor_num = 1)
#' plot_factor(sim_object = output_obj, factor_num = 'all', type = 'histogram')
#' @export
plot_factor <- function(sim_object = NULL, factor_num = NULL, type = "scatter", show.legend = TRUE) {
  if (is.null(sim_object) || is.null(factor_num)) {
    sim_object <- OmixCraftHD(vector_features = c(2000,2000), n_samples = 50, sigmas_vector = c(3,5), n_factors = 3, num.factor = 'multiple')
    factor_num <- 1
    message("Note: No input provided. Using default simulated data.\nFactor_num set to 1.")
  }

  combined_list <- c(sim_object$list_alphas, sim_object$list_gammas)
  end_digits <- sub(".*?(\\d+)$", "\\1", names(combined_list))
  unique_end_digits <- !duplicated(end_digits)
  factor_scores <- combined_list[unique_end_digits]
  names(factor_scores) <- gsub("alpha|gamma", "score", names(factor_scores))

  samples <- rownames(sim_object$concatenated_datasets[[1]])

  build_scatter <- function(df, title) {
    p <- ggplot(df, aes(x = Index, y = factor, color = factor)) +
      geom_point() +
      labs(x = "Samples", y = "Factor score", color = "") +
      scale_color_viridis_c() +
      theme_bw() +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
    if (!show.legend) p <- p + theme(legend.position = "none")
    return(p)
  }

  build_histogram <- function(df, title) {
    p <- ggplot(df, aes(x = factor, fill = after_stat(x))) +
      geom_histogram(bins = 100) +
      labs(x = "Scores", y = "Frequency", fill = "") +
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

  if (factor_num == "all") {
    plot_list <- list()
    for (i in seq_along(factor_scores)) {
      score <- factor_scores[[i]]
      df <- data.frame(samples = samples, factor = score)
      df$Index <- as.numeric(gsub("\\D", "", df$samples))
      df <- df[order(df$Index), ]
      title <- paste("Factor", i)
      p <- if (type == "scatter") build_scatter(df, title) else build_histogram(df, title)
      plot_list[[i]] <- p
    }
    return(gridExtra::grid.arrange(grobs = plot_list, ncol = length(plot_list)))
  } else if (factor_num %in% seq_along(factor_scores)) {
    score <- factor_scores[[factor_num]]
    df <- data.frame(samples = samples, factor = score)
    df$Index <- as.numeric(gsub("\\D", "", df$samples))
    df <- df[order(df$Index), ]
    title <- paste("Factor", factor_num)
    return(if (type == "scatter") build_scatter(df, title) else build_histogram(df, title))
  } else {
    stop("Invalid factor_num. Use a number within available range or 'all'.")
  }
}

utils::globalVariables(c("Index", "x"))
