#' @name plot_simData
#' @title Visualizing the simulated data using image map and 3D visualization
#' @param sim_object R object containing simulated data to be plotted
#' @param type type of the plot. Heatmap for image plot and 3D for persp 3D plot
#' @importFrom graphics image persp
#' @examples
#' # Examples
#' output_obj <- OmixCraftHD(
#'   vector_features = c(2000,3000),
#'   sigmas_vector=c(8,5),
#'   n_samples=100,
#'   n_factors=5,
#'   num.factor='multiple',
#'   advanced_dist='mixed'
#' )
#' plot_simData(sim_object = output_obj, type = "heatmap")
#' plot_simData(sim_object = output_obj, type = "3D")
#' @export
plot_simData <- function(sim_object, type='heatmap'){
  datasets <- sim_object$concatenated_datasets
  var_sigma <- sim_object$var_sigma
  main_font_size = 1
  # Extract datasets based on iteration and var_sigma
  dataset <- datasets[[1]]

  if (type == "heatmap") {
    # Plot image
    image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples",main = '', cex.main = main_font_size)#main = paste("Noise (var.):", var_sigma), cex.main = main_font_size)
  } else if (type == "3D") {
    # Plot persp
    persp(c(1:dim(t(dataset))[2]), c(1:dim(t(dataset))[1]), dataset, theta = 325, phi = 15, col = "yellow", xlab = "", ylab = "Samples", zlab = " ")#, main = paste("Noise (var.):", var_sigma, ", iter:", iteration), cex.main = main_font_size)
  } else {
    stop("Invalid type. Choose 'heatmap' or '3D'.")
  }
}
#plot_simData(multi_factor_output, type = "heatmap")
#plot_simData(data_single_unique, type = "3D")#, var_sigma = 15, iteration = 1)

# Type of the plot can be 'Heatmap' or '3D'
# Heatmap is generated using the
