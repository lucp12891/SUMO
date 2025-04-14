#' @name plot_simData
#' @title Visualizing the simulated data using heatmap or 3D surface plot
#' @description Generates a visual representation of the simulated omics data either as a heatmap or a 3D surface plot. You can select which dataset to visualize: the merged/concatenated matrix, omic one, or omic two.
#'
#' @param sim_object R object containing simulated data as returned by `OmixCraftHD`.
#' @param data Character. Specifies which data matrix to visualize. Options are "merged" (or "concatenated"), "omic.one", or "omic.two".
#' @param type Character. Type of plot: either "heatmap" for a 2D image plot or "3D" for a 3D perspective surface plot.
#' @importFrom graphics image persp abline
#' @importFrom readxl read_excel
#' @importFrom readr read_csv
#' @importFrom stringr str_detect
#' @examples
#' output_obj <- OmixCraftHD(
#'   vector_features = c(4000,3000),
#'   n_samples=100,
#'   n_factors=2,
#'   signal.samples = NULL,
#'   signal.features.one = NULL,
#'   signal.features.two = NULL,
#'   snr = 2.5,
#'   num.factor='multiple',
#'   advanced_dist='mixed')
#'
#' plot_simData(sim_object = output_obj, data = "merged", type = "heatmap")
#' plot_simData(sim_object = output_obj, data = "omic.one", type = "3D")
#' @export
plot_simData <- function(sim_object, data = "merged", type = "heatmap") {
  main_font_size <- 1
  # length of the first dataset
  len <- ncol(sim_object$omic.one[[1]])

  # Select dataset based on input
  dataset <- switch(tolower(data),
                    "merged" = t(sim_object$concatenated_datasets[[1]]),
                    "concatenated" = t(sim_object$concatenated_datasets[[1]]),
                    "omic.one" = t(sim_object$omic.one[[1]]),
                    "omic.two" = t(sim_object$omic.two[[1]]),
                    stop("Invalid data type. Choose from 'merged', 'omic.one', or 'omic.two'.")
  )

  if (type == "heatmap") {
    image(1:ncol(dataset), 1:nrow(dataset), t(dataset[nrow(dataset):1, ]),
          xlab = "Samples", ylab = "Features", main = "", cex.main = main_font_size)

    if (tolower(data) %in% c('merged', 'concatenated')) {
      abline(h = len, col = "brown", lty = "dashed")  # Feature split line
    }

  } else if (type == "3D") {
    persp(c(1:dim(t(dataset))[2]), c(1:dim(t(dataset))[1]),
          dataset, theta = 325, phi = 15, col = "yellow",
          xlab = "", ylab = "Samples", zlab = " ")
  } else {
    stop("Invalid type. Choose 'heatmap' or '3D'.")
  }
}

