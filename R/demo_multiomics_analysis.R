#' @name demo_multiomics_analysis
#' @title Demonstration of SUMO Utility in Multi-Omics Analysis using MOFA2
#' @description Run a complete MOFA2 workflow on either SUMO-generated data or the
#'   real-world CLL dataset. The function handles preprocessing and model training
#'   (preferring MOFA2's basilisk; falling back to a user reticulate env if configured),
#'   or loads a bundled pretrained model. It then creates summary visualizations and
#'   can export a multi-slide PowerPoint report.
#'
#' @param data_type Character. "SUMO" (synthetic) or "real_world" (CLL).
#' @param export_pptx Logical. If TRUE, write a PowerPoint report (multiple slides).
#'   Default TRUE.
#' @param verbose Logical. If TRUE, print progress messages. Default TRUE.
#' @param use_pretrained One of "auto", "always", "never".
#'   * "auto": train if a backend is available, otherwise load a pretrained model.
#'   * "always": always load a pretrained model and skip training.
#'   * "never": always train (requires a working Python backend for MOFA2).
#'
#' @details
#' Backend selection. Training prefers MOFA2's basilisk backend when available;
#' otherwise a reticulate/conda environment is used if configured via
#' \code{sumo_setup_mofa()}. If neither is available and \code{use_pretrained = "auto"},
#' the function loads a pretrained model shipped under \code{inst/extdata/}.
#'
#' PowerPoint contents (when \code{export_pptx = TRUE}):
#'   - Title slide (dataset label and generation date)
#'   - Data overview: \code{plot_data_overview()}
#'   - Factor correlation: \code{plot_factor_cor()}
#'   - Variance explained:
#'       - by view x factor, and
#'       - by group x factor (with totals) via \code{plot_variance_explained()}
#'   - Factor visualizations:
#'       - beeswarms for factors 1-3 via \code{plot_factor()}
#'       - a customized F1 vs F2 plot
#'       - scatter plots of factor combinations via \code{plot_factors()}
#'   - Feature weights:
#'       - \code{plot_weights()} and \code{plot_top_weights()} (first view, factor 1)
#'   - Input-data views:
#'       - heatmap via \code{plot_data_heatmap()} and
#'       - feature-factor scatter via \code{plot_data_scatter()} (second view if present)
#'   - Non-linear embedding (if available): t-SNE via \code{run_tsne()} + \code{plot_dimred()}
#'   - Table slides (heads/summaries):
#'       - sample metadata (head)
#'       - total R^2 per view/group (head)
#'       - R^2 per factor x view/group (head)
#'       - dimensions summary (factors/weights/data)
#'       - long-format heads from \code{get_factors()}, \code{get_weights()}, \code{get_data()}
#'
#' Plots are rasterized for portability when embedding in PPT (vector export is
#' used when supported).
#'
#' @return Invisibly returns the trained (or loaded) MOFA model object.
#'
#' @examples
#' if (
#'   interactive() &&
#'   requireNamespace("MOFA2", quietly = TRUE) &&
#'   requireNamespace("systemfonts", quietly = TRUE) &&
#'   utils::packageVersion("systemfonts") >= "1.1.0" &&
#'   identical(Sys.getenv("NOT_CRAN"), "true")
#' ) {
#'   # Use pretrained models (no Python needed):
#'   demo_multiomics_analysis("SUMO",       export_pptx = TRUE, use_pretrained = "always")
#'   demo_multiomics_analysis("real_world", export_pptx = TRUE, use_pretrained = "always")
#'
#'   # To train (when basilisk or a reticulate env is available):
#'   # demo_multiomics_analysis("real_world", export_pptx = TRUE, use_pretrained = "never")
#' }
#'
#' @seealso [simulate_twoOmicsData()], [plot_factor()], [plot_weights()],
#'   [sumo_setup_mofa()], [sumo_mofa_backend()], [sumo_load_pretrained_mofa()]
#' @keywords MOFA demo multi-omics synthetic-data
#' @importFrom utils packageVersion data head capture.output
#' @export
demo_multiomics_analysis <- function(
    data_type = c("SUMO", "real_world"),
    export_pptx = TRUE,
    verbose = TRUE,
    use_pretrained = c("auto", "always", "never")
) {
  data_type     <- match.arg(data_type)
  use_pretrained <- match.arg(use_pretrained)
  log_msg <- function(msg) if (verbose) message(msg)

  # Core deps used for either path (train or load)
  required_pkgs <- c("MOFA2", "ggplot2")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package ", pkg, " is required but not installed.")
    }
  }

  # Optional PowerPoint deps
  pptx_supported <- export_pptx &&
    requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("rvg", quietly = TRUE) &&
    requireNamespace("systemfonts", quietly = TRUE) &&
    utils::packageVersion("systemfonts") >= "1.1.0"
  if (export_pptx && !pptx_supported) {
    warning("PowerPoint export requires 'officer', 'rvg', and 'systemfonts >= 1.1.0'. Disabling export.")
    export_pptx <- FALSE
  }

  # --- Backend selection (basilisk \u2192 reticulate \u2192 pretrained) ---
  backend <- sumo_mofa_backend()

  # Helper: load a pretrained model and render plots/report
  load_pretrained <- function(kind) {
    log_msg(paste0("Loading pretrained MOFA (", kind, ") from inst/extdata ..."))
    model <- sumo_load_pretrained_mofa(if (kind == "SUMO") "SUMO" else "CLL")
    generate_visuals_and_report(model, kind, export_pptx)
    invisible(model)
  }

  # Training wrapper
  run_mofa_analysis <- function(data_list, n_factors = 2, feature_names = NULL, use_basilisk = TRUE) {
    mofa_data <- if (!is.null(feature_names)) {
      MOFA2::create_mofa(data_list, feature_names = feature_names)
    } else {
      MOFA2::create_mofa(data_list)
    }

    data_opts <- MOFA2::get_default_data_options(mofa_data)
    data_opts$scale_views   <- FALSE
    data_opts$scale_groups  <- TRUE
    data_opts$center_groups <- TRUE

    model_opts <- MOFA2::get_default_model_options(mofa_data)
    model_opts$num_factors <- n_factors

    train_opts <- MOFA2::get_default_training_options(mofa_data)
    train_opts$maxiter <- 1000
    train_opts$convergence_mode <- "slow"
    train_opts$seed <- 123

    mofa_obj <- MOFA2::prepare_mofa(
      object           = mofa_data,
      data_options     = data_opts,
      model_options    = model_opts,
      training_options = train_opts
    )

    # If using reticulate and a stored env, politely activate it for this call
    if (identical(backend$method, "reticulate") &&
        !is.null(backend$envname) &&
        requireNamespace("reticulate", quietly = TRUE)) {
      try(reticulate::use_condaenv(backend$envname, required = FALSE), silent = TRUE)
    }

    outfile <- tempfile(fileext = ".hdf5")
    suppressWarnings(MOFA2::run_mofa(mofa_obj, outfile, use_basilisk = use_basilisk))
  }

  # Minimal reporting/plots (+ optional PPTX cover slide)
  # Minimal reporting/plots (+ full PPT export)
  generate_visuals_and_report <- function(model, data_label, export_pptx, output_name = NULL) {
    # ---- Ensure some demo metadata so downstream plots work ----
    Nsamples <- sum(model@dimensions$N)
    md <- try(MOFA2::samples_metadata(model), silent = TRUE)
    need_meta <- inherits(md, "try-error") || !is.data.frame(md) ||
      !"sample" %in% names(md) || !"condition" %in% names(md) || !"age" %in% names(md)
    if (need_meta) {
      sample_metadata <- data.frame(
        sample    = MOFA2::samples_names(model)[[1]],
        condition = sample(c("A","B"), size = Nsamples, replace = TRUE),
        age       = sample(1:100, size = Nsamples, replace = TRUE),
        stringsAsFactors = FALSE
      )
      MOFA2::samples_metadata(model) <- sample_metadata
    } else {
      sample_metadata <- md
    }

    # ---- Collect plots safely ----
    plots <- list(); titles <- character(0)
    add_plot <- function(p, ttl) {
      if (!inherits(p, "try-error") && !is.null(p)) {
        plots[[length(plots) + 1L]] <<- p
        titles[length(titles) + 1L] <<- ttl
      }
    }

    # Overview & factor correlation
    add_plot(MOFA2::plot_data_overview(model),                                  "Data overview")
    add_plot(try(MOFA2::plot_factor_cor(model), silent = TRUE),                  "Factor correlation")

    # Variance explained
    add_plot(MOFA2::plot_variance_explained(model, x = "view",  y = "factor"),   "Variance explained: by view \u00D7 factor")
    add_plot(MOFA2::plot_variance_explained(model, x = "group", y = "factor",
                                            plot_total = TRUE),                  "Variance explained: by group \u00D7 factor (with totals)")

    # Factors
    add_plot(try(MOFA2::plot_factor(model, factor = 1:3,
                                    color_by = "condition", shape_by = "condition"), silent = TRUE),
             "Single-factor beeswarm (F1\u2013F3)")

    p_custom <- try(MOFA2::plot_factor(model, factors = c(1,2),
                                       color_by = "condition", dot_size = 3,
                                       dodge = TRUE, legend = FALSE,
                                       add_violin = TRUE, violin_alpha = 0.25), silent = TRUE)
    if (!inherits(p_custom, "try-error")) {
      p_custom <- p_custom +
        ggplot2::scale_color_manual(values = c(A = "black", B = "red")) +
        ggplot2::scale_fill_manual(values  = c(A = "black", B = "red"))
    }
    add_plot(p_custom, "Custom factor plot (F1 vs F2)")

    add_plot(try(MOFA2::plot_factors(model, factors = 1:3, color_by = "condition"), silent = TRUE),
             "Factor scatter (combinations of F1\u2013F3)")

    # Views for weights/heatmaps/scatter
    wlist  <- try(MOFA2::get_weights(model, as.data.frame = FALSE), silent = TRUE)
    vnames <- if (!inherits(wlist, "try-error") && length(wlist)) names(wlist) else names(MOFA2::get_data(model))
    if (length(vnames) == 0) vnames <- c("view_1", "view_2")

    # Weights
    add_plot(try(MOFA2::plot_weights(model, view = vnames[1], factor = 1,
                                     nfeatures = 10, scale = TRUE, abs = FALSE), silent = TRUE),
             sprintf("Weights (view = %s, factor = 1)", vnames[1]))
    add_plot(try(MOFA2::plot_top_weights(model, view = vnames[1], factor = 1, nfeatures = 10), silent = TRUE),
             sprintf("Top weights (view = %s, factor = 1)", vnames[1]))

    # Heatmap & scatter on second view if present
    view2 <- if (length(vnames) >= 2) vnames[2] else vnames[1]
    add_plot(try(MOFA2::plot_data_heatmap(model, view = view2, factor = 1, features = 20,
                                          cluster_rows = TRUE, cluster_cols = FALSE,
                                          show_rownames = TRUE, show_colnames = FALSE), silent = TRUE),
             sprintf("Data heatmap (view = %s, factor = 1)", view2))
    add_plot(try(MOFA2::plot_data_scatter(model, view = view2, factor = 1, features = 5,
                                          add_lm = TRUE, color_by = "condition"), silent = TRUE),
             sprintf("Data scatter (view = %s, factor = 1)", view2))

    # t-SNE (optional)
    model_tsne <- try(MOFA2::run_tsne(model), silent = TRUE)
    if (!inherits(model_tsne, "try-error")) {
      add_plot(try(MOFA2::plot_dimred(model_tsne, method = "TSNE", color_by = "condition"), silent = TRUE),
               "t-SNE on MOFA factors")
    }

    # ---- Key tables/heads ----
    # variance explained (stored in cache)
    r2_total_head     <- try(utils::head(model@cache$variance_explained$r2_total[[1]]),     silent = TRUE)
    r2_per_factor_head<- try(utils::head(model@cache$variance_explained$r2_per_factor[[1]]),silent = TRUE)

    # dims summary
    factors_list <- try(MOFA2::get_factors(model, factors = "all"), silent = TRUE)
    weights_list <- try(MOFA2::get_weights(model, views = "all", factors = "all"), silent = TRUE)
    data_list    <- try(MOFA2::get_data(model), silent = TRUE)

    dims_tbl <- data.frame(Section = character(), Name = character(), Nrow = integer(), Ncol = integer(), stringsAsFactors = FALSE)
    if (!inherits(factors_list, "try-error")) {
      for (nm in names(factors_list)) dims_tbl <- rbind(dims_tbl, data.frame(Section="factors", Name=nm, Nrow=nrow(factors_list[[nm]]), Ncol=ncol(factors_list[[nm]])))
    }
    if (!inherits(weights_list, "try-error")) {
      for (nm in names(weights_list)) dims_tbl <- rbind(dims_tbl, data.frame(Section="weights", Name=nm, Nrow=nrow(weights_list[[nm]]), Ncol=ncol(weights_list[[nm]])))
    }
    if (!inherits(data_list, "try-error")) {
      for (grp in names(data_list)) for (nm in names(data_list[[grp]])) {
        dims_tbl <- rbind(dims_tbl, data.frame(Section = paste0("data:", grp), Name = nm,
                                               Nrow = nrow(data_list[[grp]][[nm]]), Ncol = ncol(data_list[[grp]][[nm]])))
      }
    }

    # long-format heads
    factors_df_head <- try(utils::head(MOFA2::get_factors(model, as.data.frame = TRUE), 3), silent = TRUE)
    weights_df_head <- try(utils::head(MOFA2::get_weights(model, as.data.frame = TRUE), 3), silent = TRUE)
    data_df_head    <- try(utils::head(MOFA2::get_data(model, as.data.frame = TRUE), 3),    silent = TRUE)

    # ---- If not exporting, print & exit ----
    if (!export_pptx) {
      for (i in seq_along(plots)) { print(plots[[i]]) }
      if (!inherits(sample_metadata, "try-error")) print(utils::head(sample_metadata, 3))
      if (!inherits(r2_total_head, "try-error"))      print(r2_total_head)
      if (!inherits(r2_per_factor_head, "try-error")) print(r2_per_factor_head)
      if (nrow(dims_tbl)) print(dims_tbl)
      if (!inherits(factors_df_head, "try-error")) print(factors_df_head)
      if (!inherits(weights_df_head, "try-error")) print(weights_df_head)
      if (!inherits(data_df_head, "try-error"))    print(data_df_head)
      return(invisible(NULL))
    }

    # ---- Build PPTX ----
    ppt <- officer::read_pptx()
    ppt <- officer::add_slide(ppt, layout = "Title Slide", master = "Office Theme")
    ppt <- officer::ph_with(ppt, value = paste("MOFA Results -", data_label),
                            location = officer::ph_location_type(type = "ctrTitle"))
    ppt <- officer::ph_with(ppt, value = paste("Generated:", Sys.Date()),
                            location = officer::ph_location_type(type = "subTitle"))

    # Helper: add plot slide by rasterizing (robust across plot types)
    add_plot_slide <- function(ppt, plot_obj, title = NULL, width_px = 1920, height_px = 1080, res = 200) {
      f <- tempfile(fileext = ".png")
      if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png(f, width = width_px, height = height_px, res = res)
      } else {
        grDevices::png(f, width = width_px, height = height_px, res = res)
      }
      try(print(plot_obj), silent = TRUE)
      grDevices::dev.off()
      ppt <- officer::add_slide(ppt, layout = "Blank", master = "Office Theme")
      if (!is.null(title)) {
        # add a small title box
        ppt <- officer::ph_with(
          ppt, value = title,
          location = officer::ph_location(left = 0.4, top = 0.2, width = 12, height = 0.6)
        )
      }
      ppt <- officer::ph_with(ppt, officer::external_img(f), location = officer::ph_location_fullsize())
      ppt
    }

    # Helper: add a small table slide
    add_table_slide <- function(ppt, df, title) {
      ppt <- officer::add_slide(ppt, layout = "Title and Content", master = "Office Theme")
      ppt <- officer::ph_with(ppt, value = title, location = officer::ph_location_type(type = "title"))
      if (requireNamespace("flextable", quietly = TRUE)) {
        ft <- flextable::flextable(df); ft <- flextable::autofit(ft)
        ppt <- officer::ph_with(ppt, value = ft, location = officer::ph_location_type(type = "body"))
      } else {
        txt <- paste(capture.output(print(utils::head(df, 20))), collapse = "\n")
        ppt <- officer::ph_with(ppt, value = txt, location = officer::ph_location_type(type = "body"))
      }
      ppt
    }

    # Plot slides
    for (i in seq_along(plots)) {
      ppt <- add_plot_slide(ppt, plots[[i]], title = titles[i])
    }

    # Table slides
    if (!inherits(sample_metadata, "try-error"))
      ppt <- add_table_slide(ppt, utils::head(sample_metadata, 3), "Sample metadata (head)")
    if (!inherits(r2_total_head, "try-error"))
      ppt <- add_table_slide(ppt, as.data.frame(r2_total_head), "R\u00B2 total per view/group (head)")
    if (!inherits(r2_per_factor_head, "try-error"))
      ppt <- add_table_slide(ppt, as.data.frame(r2_per_factor_head), "R\u00B2 per factor \u00D7 view/group (head)")
    if (nrow(dims_tbl))
      ppt <- add_table_slide(ppt, dims_tbl, "Dimensions summary (factors / weights / data)")
    if (!inherits(factors_df_head, "try-error"))
      ppt <- add_table_slide(ppt, factors_df_head, "Factors (long format, head)")
    if (!inherits(weights_df_head, "try-error"))
      ppt <- add_table_slide(ppt, weights_df_head, "Weights (long format, head)")
    if (!inherits(data_df_head, "try-error"))
      ppt <- add_table_slide(ppt, data_df_head, "Data (long format, head)")

    # Save PPT
    outfile <- ifelse(is.null(output_name), paste0("MOFA_", data_label, ".pptx"), output_name)
    print(ppt, target = outfile)
    message("PowerPoint saved as ", normalizePath(outfile, winslash = "/", mustWork = FALSE))
    invisible(NULL)
  }


  # ===== Strategy selection =====
  if (use_pretrained == "always") {
    return(load_pretrained(if (data_type == "SUMO") "SUMO" else "CLL"))
  }
  if (use_pretrained == "auto" && identical(backend$method, "pretrained")) {
    return(load_pretrained(if (data_type == "SUMO") "SUMO" else "CLL"))
  }

  # ===== Data branches =====
  if (data_type == "SUMO") {
    log_msg("Processing SUMO-generated data...")

    set.seed(934)
    sim <- SUMO::simulate_twoOmicsData(
      vector_features = c(4000, 3000),
      n_samples       = 100,
      snr             = 0.5,
      multiomics_demo = "SHOW"
    )

    omic_data <- list(
      first_omic  = t(as.matrix(sim$concatenated_datasets[, (sim$n_features_two + 1):(sim$n_features_two + sim$n_features_one)])),
      second_omic = t(as.matrix(sim$concatenated_datasets[, 1:sim$n_features_two]))
    )
    feat_names <- list(colnames(omic_data[[1]]), colnames(omic_data[[2]]))

    model <- run_mofa_analysis(
      omic_data,
      feature_names = feat_names,
      use_basilisk  = isTRUE(backend$use_basilisk)
    )

    # Add sample metadata indicating simulated signal membership
    sample_meta <- data.frame(sample = paste0("sample_", seq_len(sim$n_samples)))
    for (i in seq_along(sim$indices_samples)) {
      sample_meta[[paste0("signal_Factor", i)]] <- seq_len(sim$n_samples) %in% sim$indices_samples[[i]]
    }
    MOFA2::samples_metadata(model) <- sample_meta

    generate_visuals_and_report(model, "SUMO", export_pptx)
    invisible(model)

  } else {
    log_msg("Processing real-world CLL data...")

    if (!requireNamespace("MOFAdata", quietly = TRUE))
      stop("Package MOFAdata is required for 'real_world' demo.")

    utils::data("CLL_data", package = "MOFAdata", envir = environment())

    # Keep 2 views for speed; coerce, clean, and ensure dimnames
    view_ids  <- c(2, 3)
    CLL_small <- get("CLL_data", envir = environment())[view_ids]

    cll <- lapply(CLL_small, function(x) {
      x <- as.matrix(x)
      x[is.na(x)] <- 0
      if (is.null(rownames(x))) rownames(x) <- paste0("F", seq_len(nrow(x)))
      if (is.null(colnames(x))) colnames(x) <- paste0("S", seq_len(ncol(x)))
      x
    })
    names(cll) <- names(CLL_small)

    stopifnot(length(cll) >= 1L)
    stopifnot(all(vapply(cll, nrow, 1L) > 0),
              all(vapply(cll, ncol, 1L) > 0))

    model <- run_mofa_analysis(
      cll,
      use_basilisk = isTRUE(backend$use_basilisk)
    )

    generate_visuals_and_report(model, "real_world", export_pptx)
    invisible(model)
  }
}
