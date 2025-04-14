#' @name demo_multiomics_analysis
#' @title Demonstration of SUMO Utility in Multi-Omics Analysis using MOFA2
#' @description Run a complete MOFA2-based analysis pipeline using either SUMO-generated or real-world CLL multi-omics data.
#' This function includes preprocessing, MOFA model training, variance decomposition visualization, and optional PowerPoint report generation.
#'
#' @param data_type Character. Options are `"SUMO"` for synthetic data or `"real_world"` for the CLL dataset.
#' @param export_pptx Logical. If `TRUE`, saves a PowerPoint summary of the analysis. Default is `TRUE`.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default is `TRUE`.
#'
#' @details The function checks for required packages such as MOFA2, MOFAdata, officer, basilisk, and others using `requireNamespace()`.
#' All downstream analysis is conditionally executed based on package availability.
#'
#' @return Invisibly returns the trained MOFA model object. Optionally saves visualizations and a PowerPoint report to disk.
#'
#' @examples
#' if (
#'   requireNamespace("MOFA2", quietly = TRUE) &&
#'   requireNamespace("MOFAdata", quietly = TRUE) &&
#'   identical(Sys.getenv("NOT_CRAN"), "true")
#' ) {
#'   demo_multiomics_analysis("SUMO", export_pptx = FALSE)
#'   demo_multiomics_analysis("real_world", export_pptx = FALSE)
#' }
#' @seealso [OmixCraftHD()], [plot_factor()], [plot_weights()]
#' @keywords multi-omics MOFA demo synthetic-data
#' @export
demo_multiomics_analysis <- function(data_type = c("SUMO", "real_world"), export_pptx = TRUE, verbose = TRUE) {
  data_type <- match.arg(data_type)

  required_pkgs <- c("readr", "officer", "readxl", "dplyr", "stringr",
                     "basilisk", "MOFA2", "data.table", "ggplot2", "MOFAdata", "rvg", "tidyverse", "grid")

  lapply(required_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop(paste("Package", pkg, "is required."))
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  })

  # Helper: log messages
  log_msg <- function(msg) if (verbose) message(msg)

  # Prepare SUMO-generated data
  prepare_sumo_data <- function() {
    set.seed(934)
    sim_object <- OmixCraftHD(vector_features = c(4000, 3000), n_samples = 100, snr = 0.5, multiomics_demo = "SHOW")
    list(
      data = sim_object$concatenated_datasets,
      noise = sim_object$noise,
      sd_smp = sim_object$sd_smp,
      n_features = list(one = sim_object$n_features_one, two = sim_object$n_features_two)
    )
  }

  # Prepare real-world CLL data
  prepare_cll_data <- function() {
    utils::data("CLL_data", package = "MOFAdata", envir = environment())
    cll_data <- CLL_data[c(2, 3)]

    cll_data <- lapply(cll_data, function(df) {
      df <- as.data.frame(df)
      df[is.na(df)] <- 0
      as.matrix(df)
    })

    names(cll_data) <- names(CLL_data)[c(2, 3)]

    list(
      data = cll_data
    )
  }

  # Run MOFA analysis
  run_mofa_analysis <- function(data_list, n_factors = 2, feature_names = NULL, metadata_url = NULL) {
    if (!is.null(feature_names)) {
      mofa_data <- MOFA2::create_mofa(data_list, use_basilisk = TRUE, feature_names = feature_names)
    } else {
      mofa_data <- MOFA2::create_mofa(data_list)
    }

    data_opts <- MOFA2::get_default_data_options(mofa_data)
    data_opts$scale_views <- FALSE
    data_opts$scale_groups <- TRUE
    data_opts$center_groups <- TRUE

    model_opts <- MOFA2::get_default_model_options(mofa_data)
    model_opts$num_factors <- n_factors

    train_opts <- MOFA2::get_default_training_options(mofa_data)
    train_opts$maxiter <- 1000
    train_opts$convergence_mode <- "slow"
    train_opts$seed <- 123

    mofa_obj <- MOFA2::prepare_mofa(
      object = mofa_data,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )

    outfile <- tempfile(fileext = ".hdf5")
    trained_model <- suppressWarnings(MOFA2::run_mofa(mofa_obj, outfile, use_basilisk = TRUE))

    if (!is.null(metadata_url)) {
      metadata <- data.table::fread(metadata_url)
      MOFA2::samples_metadata(trained_model) <- metadata
    }

    return(trained_model)
  }

  # Plot and report generation
  generate_visuals_and_report <- function(model, data_label, export_pptx, output_name = NULL) {
    MOFA2::plot_variance_explained(model, x = "group", y = "factor", plot_total = TRUE)
    MOFA2::plot_factor_cor(model)

    if (export_pptx) {
      ppt <- officer::read_pptx()
      ppt <- officer::add_slide(ppt, layout = "Title Slide", master = "Office Theme")
      ppt <- officer::ph_with(ppt, value = paste("MOFA Results -", data_label),
                              location = officer::ph_location_type(type = "ctrTitle"))
      ppt <- officer::ph_with(ppt, value = paste("Generated:", Sys.Date()),
                              location = officer::ph_location_type(type = "subTitle"))
      outfile <- ifelse(is.null(output_name), paste0("MOFA_", data_label, ".pptx"), output_name)
      print(ppt, target = outfile)
      message(paste("PowerPoint saved as", outfile))
    }
  }

  # Main pipeline
  if (data_type == "SUMO") {
    log_msg("Processing SUMO-generated data...")
    sumo <- prepare_sumo_data()
    omic_data <- list(
      first_omic = as.matrix(t(sumo$data[, (sumo$n_features$two + 1):(sumo$n_features$two + sumo$n_features$one)])),
      second_omic = as.matrix(t(sumo$data[, 1:sumo$n_features$two]))
    )
    feature_names <- list(
      colnames(t(sumo$data[, (sumo$n_features$two + 1):(sumo$n_features$two + sumo$n_features$one)])),
      colnames(t(sumo$data[, 1:sumo$n_features$two]))
    )
    model <- run_mofa_analysis(omic_data, feature_names = feature_names)
    generate_visuals_and_report(model, "SUMO", export_pptx)
    invisible(model)

  } else if (data_type == "real_world") {
    log_msg("Processing real-world CLL data...")
    cll <- prepare_cll_data()
    model <- run_mofa_analysis(cll$data)
    generate_visuals_and_report(model, "real_world", export_pptx)
    invisible(model)
  }
}
