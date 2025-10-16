#' Internal infix: "use left-hand side unless NULL, otherwise right-hand side"
#'
#' Lightweight helper to express a common defaulting pattern:
#' `a %||% b` returns `a` if it is not `NULL`, otherwise returns `b`.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Internal predicate: does an object already match the multi-omics schema?
#'
#' Checks whether a list-like object appears to follow the current
#' `simulateMultiOmics`-style schema (i.e., contains `omics`, `list_betas`,
#' and a `signal_annotation` entry with nested `features`).
#'
#' @param x A list-like object to inspect.
#' @return `TRUE` if the object looks like a standardized multi-omics object,
#'   otherwise `FALSE`.
#'
#' @keywords internal
#' @noRd
.is_multiomics_like <- function(x) {
  is.list(x) && !is.null(x$omics) && !is.null(x$list_betas) &&
    !is.null(x$signal_annotation) && !is.null(x$signal_annotation$features)
}


#' Coerce legacy simulation outputs to the current multi-omics schema
#' @name as_multiomics
#' @title Convert legacy objects (e.g., from `simulate_twoOmicsData()`) to the
#' current standardized structure used by downstream tools.
#'
#' @description
#' Normalizes outputs that may contain fields like `omic.one`, `omic.two`,
#' `list_betas`/`beta`, and `list_deltas`/`delta` into a unified structure with
#' `omics`, `list_betas` (per-omic), `signal_annotation` (samples and features),
#' and a `factor_map`. If the input already matches the current schema, it is
#' returned unchanged.
#'
#' @param x A list-like legacy simulation object (e.g., produced by
#'   `simulate_twoOmicsData()` or similar helpers). May also be a list with
#'   element `omics` that is itself a list of matrices.
#'
#' @return A standardized list with components:
#' \itemize{
#'   \item `concatenated_datasets` — list of one matrix `cbind(omic1, omic2, ...)`.
#'   \item `omics` — named list of matrices (samples x features).
#'   \item `list_alphas` — per-factor sample scores (named `alpha1`, `alpha2`, ...).
#'   \item `list_betas` — list per-omic with per-factor loadings (named `beta1`, `beta2`, ...).
#'   \item `signal_annotation` — list with `samples` and `features` indices per factor and per-omic.
#'   \item `factor_structure` — best-effort label: "shared", "unique", or "mixed".
#'   \item `factor_map` — which omics each factor loads on.
#' }
#'
#' @export
as_multiomics <- function(x) {
  # local helpers (avoid relying on %||%)
  first_non_null <- function(...) {
    for (y in list(...)) if (!is.null(y)) return(y)
    NULL
  }
  is_mat_like <- function(m) is.matrix(m) || inherits(m, "Matrix")

  # --- 1) Gather omics matrices (modern or legacy) ---
  omics <- NULL

  # Modern: x$omics is a list of matrices
  if (is.list(x$omics) && length(x$omics) > 0) {
    omics <- x$omics
  } else {
    # Legacy: simulate_twoOmicsData() style
    omic1 <- x$omic.one
    omic2 <- x$omic.two
    if (is.list(omic1)) omic1 <- omic1[[1]]
    if (is.list(omic2)) omic2 <- omic2[[1]]

    # Fallback: split concatenated matrix if available
    if (is.null(omic1) || is.null(omic2)) {
      cd <- x$concatenated_datasets
      if (is.list(cd)) cd <- cd[[1]]
      if (!is.null(cd) && is.matrix(cd)) {
        n1 <- first_non_null(x$n_features_one,
                             if (!is.null(colnames(cd))) sum(grepl("^omic1_", colnames(cd))) else NULL)
        n2 <- first_non_null(x$n_features_two,
                             if (!is.null(colnames(cd))) sum(grepl("^omic2_", colnames(cd))) else NULL)
        if (!is.null(n1) && !is.null(n2) && (n1 + n2) == ncol(cd)) {
          left_is_omic1 <- if (!is.null(colnames(cd))) grepl("^omic1_", colnames(cd))[1] else TRUE
          if (!isTRUE(left_is_omic1)) {
            omic2 <- cd[, seq_len(n2), drop = FALSE]
            omic1 <- cd[, (n2 + 1):(n2 + n1), drop = FALSE]
          } else {
            omic1 <- cd[, seq_len(n1), drop = FALSE]
            omic2 <- cd[, (n1 + 1):(n1 + n2), drop = FALSE]
          }
        }
      }
    }

    omics <- list(omic1 = omic1, omic2 = omic2)
  }

  # Keep only non-NULL entries
  if (is.null(omics)) omics <- list()
  omics <- omics[!vapply(omics, is.null, logical(1))]

  # Also support the case x itself is a list of matrices
  if (length(omics) == 0L && is.list(x) && length(x) > 0 &&
      all(vapply(x, is_mat_like, logical(1)))) {
    omics <- x
  }

  if (length(omics) == 0L)
    stop("as_multiomics(): could not locate omics matrices in `x` (expected x$omics, x$omic.one/x$omic.two, or x$concatenated_datasets).")

  # Coerce to numeric matrices, basic checks, and ensure names
  for (i in seq_along(omics)) {
    m <- as.matrix(omics[[i]])
    storage.mode(m) <- "double"
    if (length(dim(m)) != 2L) stop("Each omic must be a 2D matrix.")
    # Assign names if missing
    if (is.null(rownames(m))) rownames(m) <- paste0("sample_", seq_len(nrow(m)))
    if (is.null(colnames(m))) {
      nm <- names(omics)[i]
      if (is.null(nm) || !nzchar(nm)) nm <- paste0("omic", i)
      colnames(m) <- paste0(nm, "_feature_", seq_len(ncol(m)))
    }
    omics[[i]] <- m
  }

  # Ensure consistent sample counts
  nrows <- vapply(omics, nrow, integer(1))
  if (length(unique(nrows)) != 1L)
    stop("All omics must have the same number of rows (samples).")

  # Set names if missing
  if (is.null(names(omics)) || any(!nzchar(names(omics))))
    names(omics) <- paste0("omic", seq_along(omics))

  # --- 2) list_alphas (sample-level factors) ---
  list_alphas <- first_non_null(x$list_alphas, x$alpha)
  if (!is.null(list_alphas)) {
    if (!is.list(list_alphas)) list_alphas <- NULL
    else {
      if (is.null(names(list_alphas)) || any(!grepl("^alpha\\d+$", names(list_alphas))))
        names(list_alphas) <- paste0("alpha", seq_along(list_alphas))
    }
  }

  # --- 3) list_betas: per-omic feature loadings per factor ---
  # Accept either a per-omic list-of-lists or legacy betas/deltas split
  if (is.list(x$list_betas) && all(vapply(x$list_betas, is.list, logical(1))) && !is.null(names(x$list_betas))) {
    list_betas <- x$list_betas
  } else {
    betas1 <- first_non_null(x$list_betas, x$beta)
    betas2 <- first_non_null(x$list_deltas, x$delta)
    if (!is.null(betas1)) {
      if (is.null(names(betas1))) names(betas1) <- paste0("beta", seq_along(betas1))
      else names(betas1) <- sub("^delta(\\d+)$", "beta\\1", sub("^beta(\\d+)$", "beta\\1", names(betas1)))
    }
    if (!is.null(betas2)) {
      if (is.null(names(betas2))) names(betas2) <- paste0("beta", seq_along(betas2))
      else names(betas2) <- sub("^delta(\\d+)$", "beta\\1", sub("^beta(\\d+)$", "beta\\1", names(betas2)))
    }
    list_betas <- list(omic1 = betas1, omic2 = betas2)
  }
  # Keep only entries that correspond to present omics
  list_betas <- list_betas[intersect(names(list_betas), names(omics))]

  # --- 4) signal_annotation (samples + feature indices per omic/factor) ---
  samples_annot <- first_non_null(
    if (!is.null(x$signal_annotation)) x$signal_annotation$samples else NULL,
    x$indices_samples, x$assigned_indices_samples
  )

  feat_annot <- setNames(vector("list", length(omics)), names(omics))
  for (onm in names(feat_annot)) feat_annot[[onm]] <- list()

  # Fill feature indices from explicit annotation or infer from nonzero betas
  for (onm in names(list_betas)) {
    blist <- list_betas[[onm]]
    if (is.null(blist)) next
    for (bn in names(blist)) {
      fid <- sub("^beta", "", bn)
      explicit_idx <- NULL
      if (!is.null(x$signal_annotation) && !is.null(x$signal_annotation[[onm]]))
        explicit_idx <- x$signal_annotation[[onm]][[bn]]
      idx <- if (!is.null(explicit_idx)) as.integer(explicit_idx) else which(abs(blist[[bn]]) > 0)
      feat_annot[[onm]][[paste0("factor", fid)]] <- as.integer(idx)
    }
  }

  signal_annotation <- list(
    samples  = samples_annot,
    features = feat_annot
  )

  # --- 5) factor_map (which omics each factor hits) ---
  factor_names <- character(0)
  for (onm in names(list_betas)) {
    nm <- names(list_betas[[onm]])
    if (!is.null(nm)) factor_names <- union(factor_names, nm)
  }
  factor_ids <- sort(unique(as.integer(gsub("\\D", "", factor_names))))
  factor_map <- setNames(vector("list", length(factor_ids)), paste0("factor", factor_ids))
  for (f in factor_ids) {
    hits <- integer(0)
    for (onm in names(list_betas)) {
      nm <- names(list_betas[[onm]])
      if (!is.null(nm) && any(grepl(paste0("^beta", f, "$"), nm))) {
        hits <- c(hits, which(names(omics) == onm))
      }
    }
    factor_map[[paste0("factor", f)]] <- hits
  }

  # --- 6) concatenated_datasets (as a list of one matrix) ---
  merged <- do.call(cbind, omics)
  concatenated_datasets <- list(merged)

  # --- 7) best-effort factor_structure label ---
  if (length(factor_map) == 0L) {
    factor_structure <- NA_character_
  } else {
    hit_counts <- vapply(factor_map, length, integer(1))
    all_hit <- length(unique(hit_counts)) == 1L && unique(hit_counts) == length(omics)
    if (all_hit) factor_structure <- "shared"
    else if (all(hit_counts %in% c(0L, 1L)) && any(hit_counts == 1L)) factor_structure <- "unique"
    else factor_structure <- "mixed"
  }

  list(
    concatenated_datasets = concatenated_datasets,
    omics             = omics,
    list_alphas       = list_alphas,
    list_betas        = list_betas,
    signal_annotation = signal_annotation,
    factor_structure  = factor_structure,
    factor_map        = factor_map
  )
}
