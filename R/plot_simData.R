#' @name plot_simData
#' @title Visualizing the simulated data using heatmap or 3D surface plot
#' @description Generates a visual representation of the simulated omics data either as a heatmap or a 3D surface plot.
#' You can select which dataset to visualize: the merged/concatenated matrix, or any individual omic (e.g., "omic1", "omic2", etc.).
#'
#' @param sim_object R object containing simulated data as returned by `simulate_twoOmicsData` and `simulateMultiOmics`.
#' @param data Character. Specifies which data matrix to visualize. Options are "merged" (or "concatenated"), "omic.one", or "omic.two".
#' @param type Character. Type of plot: either "heatmap" for a 2D image plot or "3D" for a 3D perspective surface plot.
#' @importFrom graphics image persp abline rect legend
#' @importFrom readxl read_excel
#' @importFrom readr read_csv
#' @importFrom stringr str_detect
#' @examples
# Simulate data
#' output_obj <- simulateMultiOmics(
#'   vector_features = c(3000, 2500, 2000),
#'   n_samples = 100,
#'   n_factors = 3,
#'   snr = 3,
#'   signal.samples = c(5, 1),
#'   signal.features = list(
#'     c(3, 0.3),
#'     c(2.5, 0.25),
#'     c(2, 0.2)
#'   ),
#'   factor_structure = "mixed",
#'   num.factor = "multiple",
#'   seed = 123
#' )
#'
#' # Visualize merged heatmap
#' plot_simData(sim_object = output_obj, data = "merged", type = "heatmap")
#'
#' # Visualize omic2 in 3D
#' plot_simData(sim_object = output_obj, data = "omic2", type = "heatmap")
#'
#' @export
plot_simData <- function(sim_object, data = "merged", type = "heatmap") {
  main_font_size <- 1
  data_lower <- tolower(data)
  dataset <- NULL
  feature_split <- NULL

  if (data_lower %in% c("merged", "concatenated")) {
    dataset <- sim_object$concatenated_datasets
    if (is.list(dataset)) dataset <- dataset[[1]]
    dataset <- t(dataset)
    if (!is.null(sim_object$omics)) {
      feature_split <- cumsum(sapply(sim_object$omics, ncol))
    }
  } else if (!is.null(sim_object$omics) && data %in% names(sim_object$omics)) {
    dataset <- sim_object$omics[[data]]
    dataset <- t(dataset)
  } else {
    stop("Invalid `data` argument.")
  }

  if (!is.matrix(dataset)) {
    dataset <- as.matrix(dataset)
  }

  if (type == "heatmap") {
    image(1:ncol(dataset), 1:nrow(dataset), t(dataset[nrow(dataset):1, ]),
          xlab = "Samples", ylab = "Features", main = "", cex.main = main_font_size)

    if (!is.null(feature_split) && data_lower %in% c("merged", "concatenated")) {
      for (h in feature_split[-length(feature_split)]) {
        abline(h = h, col = "brown", lty = "dashed")
      }
    }

    #} else if (type == "3D") {
    #  if (!is.matrix(dataset) || !is.numeric(dataset)) {
    #    stop("3D plotting requires a numeric matrix.")
    #  }
    #  if (nrow(dataset) < 2 || ncol(dataset) < 2) {
    #    stop("3D plot requires at least a 2x2 matrix.")
    #  }
    #  if (anyNA(dataset)) {
    #    stop("Dataset contains NA values which are not allowed in 3D plots.")
    #  }
    #
    #  persp(x = 1:ncol(dataset), y = 1:nrow(dataset), z = dataset,
    #        theta = 325, phi = 15, col = "yellow",
    #        xlab = "", ylab = "Samples", zlab = "")
  } else {
    stop("Invalid `type`. Choose 'heatmap'. 3D plot not support by this version now.")
  }
}
