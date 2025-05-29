#' @name demo_multiomics_analysis
#' @title Demonstration of SUMO Utility in Multi-Omics Analysis using MOFA2
#' @description Run a complete MOFA2-based analysis pipeline using either SUMO-generated or real-world CLL multi-omics data.
#' This function includes preprocessing, MOFA model training, variance decomposition visualization, and optional PowerPoint report generation.
#'
#' @param data_type Character. Options are `"SUMO"` for synthetic data or `"real_world"` for the CLL dataset.
#' @param export_pptx Logical. If `TRUE`, saves a PowerPoint summary of the analysis. Default is `TRUE`.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default is `TRUE`.
#'
#' @details PowerPoint generation is skipped if required packages (`officer`, `rvg`, and `systemfonts >= 1.1.0`) are not available.
#'
#' @return Invisibly returns the trained MOFA model object.
#'
#' @examples
#' if (
#'   requireNamespace("MOFA2", quietly = TRUE) &&
#'   requireNamespace("MOFAdata", quietly = TRUE) &&
#'   requireNamespace("systemfonts", quietly = TRUE) &&
#'   utils::packageVersion("systemfonts") >= "1.1.0" &&
#'   identical(Sys.getenv("NOT_CRAN"), "true")
#' ) {
#'   demo_multiomics_analysis("SUMO", export_pptx = FALSE)
#'   demo_multiomics_analysis("real_world", export_pptx = FALSE)
#' }
#'
#' @seealso [simulate_twoOmicsData()], [plot_factor()], [plot_weights()]
#' @keywords MOFA demo multi-omics synthetic-data
#' @importFrom utils packageVersion data
#' @importFrom data.table fread
#' @export
demo_multiomics_analysis <- function(data_type = c("SUMO", "real_world"), export_pptx = TRUE, verbose = TRUE) {
  data_type <- match.arg(data_type)

  # Check core dependencies
  required_pkgs <- c("MOFA2", "MOFAdata", "data.table", "ggplot2", "basilisk")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  # Check for PowerPoint support
  pptx_supported <- export_pptx &&
    requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("rvg", quietly = TRUE) &&
    requireNamespace("systemfonts", quietly = TRUE) &&
    utils::packageVersion("systemfonts") >= "1.1.0"

  if (export_pptx && !pptx_supported) {
    warning("PowerPoint export requires 'officer', 'rvg', and 'systemfonts >= 1.1.0'. Disabling export.")
    export_pptx <- FALSE
  }

  log_msg <- function(msg) if (verbose) message(msg)

  prepare_sumo_data <- function() {
    set.seed(934)
    sim_object <- simulate_twoOmicsData(
      vector_features = c(4000, 3000),
      n_samples = 100,
      snr = 0.5,
      multiomics_demo = "SHOW"
    )
    list(
      data = sim_object$concatenated_datasets,
      n_features = list(one = sim_object$n_features_one, two = sim_object$n_features_two)
    )
  }

  prepare_cll_data <- function() {
    utils::data("CLL_data", package = "MOFAdata", envir = environment())
    CLL_data <- get("CLL_data", envir = environment())
    cll_data <- CLL_data[c(2, 3)]

    cll_data <- lapply(cll_data, function(df) {
      df <- as.data.frame(df)
      df[is.na(df)] <- 0
      as.matrix(df)
    })

    names(cll_data) <- names(CLL_data)[c(2, 3)]

    list(data = cll_data)
  }

  run_mofa_analysis <- function(data_list, n_factors = 2, feature_names = NULL) {
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
    return(trained_model)
  }

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

  # === MAIN LOGIC ===
  if (data_type == "SUMO") {
    log_msg("Processing SUMO-generated data...")
    sumo <- prepare_sumo_data()
    omic_data <- list(
      first_omic = as.matrix(t(sumo$data[, (sumo$n_features$two + 1):(sumo$n_features$two + sumo$n_features$one)])),
      second_omic = as.matrix(t(sumo$data[, 1:sumo$n_features$two]))
    )
    feature_names <- list(
      colnames(omic_data[[1]]),
      colnames(omic_data[[2]])
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
