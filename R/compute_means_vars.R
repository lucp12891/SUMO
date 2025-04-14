#' @name compute_means_vars
#' @title Compute Summary Statistics for a List of Datasets
#' @description Computes overall, row-wise, and column-wise means and standard deviations for each dataset in a list. Also provides average statistics across datasets.
#'
#' @param data_list A list of numeric matrices or data frames. Each entry should be a matrix or data frame with numeric values.
#'
#' @return A named list containing:
#' \itemize{
#'   \item Overall mean and SD for each dataset.
#'   \item Average row-wise mean and SD.
#'   \item Average column-wise mean and SD.
#'   \item `mean_smp`: Average row-wise mean across all datasets.
#'   \item `sd_smp`: Average row-wise SD across all datasets.
#' }
#'
#' @importFrom dplyr mutate across
#' @importFrom stats sd
#' @importFrom magrittr %>%
#' @examples
#' # Example using simulated matrices
#' set.seed(123)
#' dataset1 <- matrix(rnorm(100, mean = 5, sd = 2), nrow = 10, ncol = 10)
#' dataset2 <- matrix(rnorm(100, mean = 10, sd = 3), nrow = 10, ncol = 10)
#' data_list <- list(dataset1, dataset2)
#' results <- compute_means_vars(data_list)
#' print(results)
#'
#' \dontrun{
#' # Example using real experimental data (requires MOFAdata)
#' if (requireNamespace("MOFAdata", quietly = TRUE)) {
#'   utils::data("CLL_data", package = "MOFAdata")
#'   CLL_data2 <- CLL_data[c(2, 3)]
#'   results <- compute_means_vars(CLL_data2)
#'   print(results)
#' }
#' }
#' @export
compute_means_vars <- function(data_list) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }

  if (!is.list(data_list) || length(data_list) < 2) {
    stop("Input must be a list of at least two data frames or matrices.")
  }

  suffixes <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")
  if (length(data_list) > length(suffixes)) {
    stop("Currently supports up to ten datasets.")
  }

  results <- list()

  for (i in seq_along(data_list)) {
    dataset <- data_list[[i]]
    suffix <- suffixes[i]

    cleaned_data <- as.data.frame(dataset) %>%
      dplyr::mutate(dplyr::across(everything(), ~ replace(., is.na(.), 0)))

    overall_mean <- mean(t(cleaned_data), na.rm = TRUE)
    overall_sd <- sd(t(cleaned_data), na.rm = TRUE)

    row_means <- rowMeans(t(cleaned_data), na.rm = TRUE)
    row_means_avg <- mean(row_means)
    row_sd <- apply(t(cleaned_data), 1, sd, na.rm = TRUE)
    row_sd_avg <- mean(row_sd)

    col_means <- colMeans(t(cleaned_data), na.rm = TRUE)
    col_means_avg <- mean(col_means)
    col_sd <- apply(t(cleaned_data), 2, sd, na.rm = TRUE)
    col_sd_avg <- mean(col_sd)

    results[[paste0("overall_mean.", suffix)]] <- overall_mean
    results[[paste0("overall_sd.", suffix)]] <- overall_sd
    results[[paste0("row_means.", suffix)]] <- row_means_avg
    results[[paste0("row_sd.", suffix)]] <- row_sd_avg
    results[[paste0("col_means.", suffix)]] <- col_means_avg
    results[[paste0("col_sd.", suffix)]] <- col_sd_avg
  }

  results[["mean_smp"]] <- mean(unlist(results[paste0("row_means.", suffixes[seq_along(data_list)])]))
  results[["sd_smp"]] <- mean(unlist(results[paste0("row_sd.", suffixes[seq_along(data_list)])]))

  return(results)
}
