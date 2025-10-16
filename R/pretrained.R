#' List pretrained MOFA models included with SUMO
#' @return Character vector of file names present in inst/extdata (empty if none).
#' @export
sumo_pretrained_mofa_available <- function() {
  base <- system.file("extdata", package = "SUMO")
  if (!nzchar(base)) return(character(0))
  files <- list.files(base, pattern = "\\.h5$|\\.hdf5$", ignore.case = TRUE)
  # Return only our expected names if present; otherwise return whatever .h5 files exist
  expected <- c("mofa_SAM_sumo_v1.h5", "mofa_SAM_cll_v1.h5")
  present_expected <- intersect(expected, files)
  if (length(present_expected)) present_expected else files
}

#' Path to a pretrained MOFA model shipped with SUMO
#' @param which One of "SUMO" or "CLL".
#' @return Full file path to the pretrained model.
#' @export
sumo_pretrained_mofa_path <- function(which = c("SUMO", "CLL")) {
  which <- match.arg(which)
  fname <- if (which == "SUMO") "mofa_SAM_sumo_v1.h5" else "mofa_SAM_cll_v1.h5"
  path <- system.file("extdata", fname, package = "SUMO")
  if (!nzchar(path) || !file.exists(path)) {
    avail <- paste(sumo_pretrained_mofa_available(), collapse = ", ")
    stop(
      "Pretrained model not found in inst/extdata: ", fname, "\n",
      "Available files: ", if (nzchar(avail)) avail else "<none>", "\n",
      "Make sure the files are under 'inst/extdata' in the source before installing.",
      call. = FALSE
    )
  }
  path
}

#' Load a pretrained MOFA model (no Python required)
#' @param which One of "SUMO" or "CLL".
#' @return A MOFA object loaded from the shipped HDF5 file.
#' @export
sumo_load_pretrained_mofa <- function(which = c("SUMO", "CLL")) {
  which <- match.arg(which)
  if (!requireNamespace("MOFA2", quietly = TRUE)) {
    stop("The R package 'MOFA2' is required to load pretrained models (install via Bioconductor).",
         call. = FALSE)
  }
  MOFA2::load_model(sumo_pretrained_mofa_path(which))
}
