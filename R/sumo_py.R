#' Detect and configure the MOFA2 backend for SUMO
#' @return A list with elements:
#'   \itemize{
#'     \item \code{method}: one of \code{"basilisk"}, \code{"reticulate"}, or \code{"pretrained"}.
#'     \item \code{use_basilisk}: logical; \code{TRUE} if using basilisk.
#'     \item \code{envname}: \code{NULL} or the name of a \code{conda} environment.
#'   }
#' @keywords internal
sumo_mofa_backend <- function() {
  # 1) Prefer basilisk via MOFA2 if available
  if (requireNamespace("MOFA2", quietly = TRUE) &&
      requireNamespace("basilisk", quietly = TRUE)) {
    return(list(method = "basilisk", use_basilisk = TRUE, envname = NULL))
  }

  # 2) Otherwise try a configured reticulate env
  if (requireNamespace("reticulate", quietly = TRUE)) {
    envname <- sumo_get_config("mofa_condaenv")
    if (isTRUE(nzchar(envname))) {
      # only select it; don't fail if missing
      try(reticulate::use_condaenv(envname, required = FALSE), silent = TRUE)
      if (isTRUE(try(reticulate::py_module_available("mofapy2"), silent = TRUE))) {
        return(list(method = "reticulate", use_basilisk = FALSE, envname = envname))
      }
    }
    # also accept already-active Python with mofapy2
    if (isTRUE(try(reticulate::py_module_available("mofapy2"), silent = TRUE))) {
      return(list(method = "reticulate", use_basilisk = FALSE, envname = NULL))
    }
  }

  # 3) Fall back to pretrained
  list(method = "pretrained", use_basilisk = FALSE, envname = NULL)
}

#' Get/set SUMO per-user configuration
#' @keywords internal
sumo_config_path <- function() {
  dir <- tools::R_user_dir("SUMO", which = "config")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, "config.json")
}

#' @keywords internal
sumo_get_config <- function(key, default = NULL) {
  path <- sumo_config_path()
  if (!file.exists(path)) return(default)
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) "")
  if (!length(txt)) return(default)
  cfg <- tryCatch(jsonlite::fromJSON(paste(txt, collapse = "\n")), error = function(e) NULL)
  if (is.null(cfg) || is.null(cfg[[key]])) default else cfg[[key]]
}

#' @keywords internal
sumo_set_config <- function(key, value) {
  path <- sumo_config_path()
  cfg <- list()
  if (file.exists(path)) {
    old <- tryCatch(jsonlite::fromJSON(paste(readLines(path, warn = FALSE), collapse = "\n")),
                    error = function(e) NULL)
    if (is.list(old)) cfg <- old
  }
  cfg[[key]] <- value
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE), con = path, useBytes = TRUE)
  invisible(TRUE)
}

#' Interactive setup for Python 'mofapy2' via reticulate (fallback when basilisk is unavailable)
#' @param envname Name of the \code{conda} environment to create/use.
#' @param py_version Python version (e.g., "3.10").
#' @return \code{TRUE} on success (and persists the env name in SUMO user config).
#' @export
sumo_setup_mofa <- function(envname = "r-mofa2", py_version = "3.10") {
  if (!interactive()) {
    message("sumo_setup_mofa() is interactive-only. See ?sumo_setup_mofa for manual steps.")
    return(invisible(FALSE))
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Please install.packages('reticulate') first.", call. = FALSE)
  }

  # Miniconda if needed
  try(reticulate::install_miniconda(), silent = TRUE)

  # Create env if missing
  envs <- try(reticulate::conda_list()$name, silent = TRUE)
  if (!isTRUE(envname %in% envs)) {
    reticulate::conda_create(envname, packages = paste0("python=", py_version))
  }

  # Activate and install deps
  reticulate::use_condaenv(envname, required = TRUE)
  reticulate::py_install(c("mofapy2", "numpy", "scipy", "pandas", "h5py", "scikit-learn"),
                         envname = envname, pip = TRUE)

  # Persist selection for future runs
  sumo_set_config("mofa_condaenv", envname)
  message("Configured reticulate env '", envname, "'. Future SUMO runs will auto-use it when basilisk is unavailable.")
  invisible(TRUE)
}

# NOTE: Avoid load-time side effects. Do not auto-activate conda here.
# Call sumo_mofa_backend() inside the analysis function that needs MOFA2,
# and activate a user-selected env there if necessary.
