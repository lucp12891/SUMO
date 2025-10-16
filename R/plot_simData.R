#' @name plot_simData
#' @title Visualize simulated multi-omics data as a heatmap
#' @description
#' Quick visualization of simulated omics data as a base R heatmap.
#' You can plot the merged/concatenated matrix across all omics or a single
#' omic layer. Optionally permute sample and/or feature order (with a seed)
#' to conceal block structure for sanity checks.
#'
#' @details
#' The function expects `sim_object$omics` to be a **named list** of numeric
#' matrices with the **same number of rows** (samples). For `data = "merged"`
#' (or `"concatenated"`), all omic matrices are column-bound in their current
#' order (subject to optional permutation) and plotted together.
#'
#' @param sim_object List-like simulation result. Must contain
#'   `omics`, a named list of numeric matrices (samples in rows, features in columns).
#' @param data Character. Which matrix to visualize:
#'   - `"merged"` or `"concatenated"`: column-bind all omics (default: `"merged"`).
#'   - a single omic name present in `names(sim_object$omics)`.
#' @param type Character. Plot type. Currently only `"heatmap"` is supported.
#' @param permute Logical. If `TRUE`, apply permutations according to
#'   `permute_samples` and `permute_features`. Default: `FALSE`.
#' @param permute_seed Integer or `NULL`. If not `NULL`, sets RNG seed **once**
#'   for reproducible permutations. Default: `NULL`.
#' @param permute_samples Logical. If `TRUE` and `permute = TRUE`, permute sample
#'   order (rows). Default: `TRUE`.
#' @param permute_features Logical. If `TRUE` and `permute = TRUE`, permute feature
#'   order (columns). Default: `TRUE`.
#'
#' @return Invisibly returns the numeric matrix that was plotted (after any permutations).
#'
#' @seealso \code{simulateMultiOmics}
#'
#' @examples
#' set.seed(123)
#' sim_object <- simulate_twoOmicsData(
#'   vector_features = c(4000, 3000),
#'   n_samples = 100,
#'   n_factors = 2,
#'   snr = 2.5,
#'   num.factor = "multiple",
#'   advanced_dist = "mixed"
#' )
#' output_obj = as_multiomics(sim_object)
#'
#' # Merged (concatenated) heatmap
#' plot_simData(output_obj, data = "merged", type = "heatmap")
#'
#' # Single omic with reproducible permutation
#' plot_simData(output_obj, data = "omic2", permute = TRUE, permute_seed = 123)
#'
#' @export
#' @importFrom graphics image
plot_simData <- function(sim_object,
                         data = "merged",
                         type = "heatmap",
                         permute = FALSE,
                         permute_seed = NULL,
                         permute_samples = TRUE,
                         permute_features = TRUE) {

  main_font_size <- 1
  data_lower <- tolower(data)
  M <- NULL

  if (!is.null(permute_seed)) set.seed(permute_seed)

  # helper: coerce to numeric matrix
  as_mat <- function(x) {
    if (!is.matrix(x)) x <- as.matrix(x)
    mode(x) <- "numeric"
    x
  }

  # ---- Build matrix M --------------------------------------------------------
  if (data_lower %in% c("merged", "concatenated")) {
    if (is.null(sim_object$omics) || length(sim_object$omics) == 0L) {
      stop("No omics found in 'sim_object' for merged view.")
    }
    omics <- sim_object$omics

    # ensure same number of samples across omics
    n_rows <- vapply(omics, nrow, integer(1))
    if (length(unique(n_rows)) != 1L) {
      stop("All matrices in 'sim_object$omics' must have the same number of rows (samples).")
    }
    n_samples <- n_rows[[1]]

    # common sample permutation
    row_perm <- if (isTRUE(permute && permute_samples)) sample.int(n_samples) else seq_len(n_samples)

    # apply per-omic feature permutation and common sample permutation
    omic_blocks <- lapply(omics, function(omic) {
      omic <- as_mat(omic)  # samples x features
      col_perm <- if (isTRUE(permute && permute_features)) sample.int(ncol(omic)) else seq_len(ncol(omic))
      omic[row_perm, col_perm, drop = FALSE]
    })

    # merge features (column-bind)
    M <- do.call(cbind, omic_blocks)

  } else if (!is.null(sim_object$omics) && data %in% names(sim_object$omics)) {
    M <- as_mat(sim_object$omics[[data]])  # samples x features
    if (isTRUE(permute)) {
      if (isTRUE(permute_samples))  M <- M[sample.int(nrow(M)), , drop = FALSE]
      if (isTRUE(permute_features)) M <- M[, sample.int(ncol(M)), drop = FALSE]
    }
  } else {
    stop("Invalid 'data' argument. Provide 'merged'/'concatenated' or a valid omic layer name.")
  }

  # ---- Plot ------------------------------------------------------------------
  if (!identical(type, "heatmap")) {
    stop("Invalid 'type'. Choose 'heatmap'.")
  }

  # image() expects z[i, j] at (x[i], y[j]); here x = samples (rows), y = features (cols)
  graphics::image(x = seq_len(nrow(M)),
                  y = seq_len(ncol(M)),
                  z = M,
                  xlab = "Samples",
                  ylab = "Features",
                  main = "",
                  cex.main = main_font_size)

  invisible(M)
}
