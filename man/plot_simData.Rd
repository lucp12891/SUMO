% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_simData.R
\name{plot_simData}
\alias{plot_simData}
\title{Visualizing the simulated data using image map and 3D visualization}
\usage{
plot_simData(sim_object, type = "heatmap")
}
\arguments{
\item{sim_object}{R object containing simulated data to be plotted}

\item{type}{type of the plot. Heatmap for image plot and 3D for persp 3D plot}
}
\description{
Visualizing the simulated data using image map and 3D visualization
}
\examples{
# Examples
output_obj <- OmixCraftHD(
  vector_features = c(2000,3000),
  sigmas_vector=c(8,5),
  n_samples=100,
  n_factors=5,
  num.factor='multiple',
  advanced_dist='mixed'
)
plot_simData(sim_object = output_obj, type = "heatmap")
plot_simData(sim_object = output_obj, type = "3D")
}
