#' @name simulateMultiOmics
#' @title Simulation of omics with predefined single or multiple latent factors in multi-omics
#'
#' @description
#' Simulate multiple omics (>=2) datasets with predefined sample-level latent factors and corresponding feature-level signal regions.
#' Each omic has unique signal structure, noise profile, and feature space.
#'
#' @details
#' This function generates synthetic multi-omics datasets for benchmarking integrative methods. Each omic layer has its own feature distribution and noise characteristics.
#'
#' Key properties:
#' - Sample signal blocks for each latent factor are non-overlapping and randomly spaced.
#' - Feature signal blocks per omic are also non-overlapping and assigned per factor.
#' - Noise can be modeled using either standard SNR scaling (default) or real data statistics (if `real_stats = TRUE`).
#' - Omics without assigned signal factors still receive background noise.
#'
#' @param vector_features Integer vector of number of features per omic (length k for k omics).
#' @param n_samples Total number of samples.
#' @param n_factors Number of latent factors.
#' @param snr Numeric. Signal-to-noise ratio.
#' @param signal.samples Length-2 vector (mean, sd) for sample-level signal values.
#' @param signal.features List of length-k vectors (mean, sd) for each omic's feature-level signal.
#' @param factor_structure Character. One of: "shared", "unique", "mixed", "partial", "custom".
#' @param num.factor Character. Either "multiple" (default) or "single" factor mode.
#' @param seed Integer seed for reproducibility (optional).
#' @param real_stats Logical. If TRUE, noise variance and mean are derived from `real_means_vars`.
#' @param real_means_vars Optional list of named vectors per omic: c(mean=..., var=...). Required if `real_stats = TRUE`.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{omics}: List of omic matrices.
#'   \item \code{concatenated_datasets}: Merged matrix of all omics.
#'   \item \code{list_alphas}: Sample-level latent factor values.
#'   \item \code{list_betas}: Feature-level loadings for each omic and factor.
#'   \item \code{signal_annotation}: List of feature and sample signal blocks.
#'   \item \code{factor_structure}: Input parameter.
#'   \item \code{factor_map}: Map of which omics each factor affects.
#' }
#'
#' @examples
#' # Example 1: Use standard SNR scaling (default)
#' sim1 <- simulateMultiOmics(
#'   vector_features = c(3000, 2500, 2000),
#'   n_samples = 100,
#'   n_factors = 3,
#'   snr = 3,
#'   signal.samples = c(5, 1),
#'   signal.features = list(
#'     c(3, 0.05),
#'     c(2.5, 0.05),
#'     c(2, 0.05)
#'   ),
#'   factor_structure = "mixed",
#'   num.factor = "multiple",
#'   seed = 123
#' )
#' plot_simData(sim_object = sim1, data = "merged", type = "heatmap")
#'
#' # Example 2: Use real stats for noise modeling
#' sim2 <- simulateMultiOmics(
#'   vector_features = c(3000, 2500, 2000),
#'   n_samples = 100,
#'   n_factors = 3,
#'   snr = 3,
#'   signal.samples = c(5, 1),
#'   signal.features = list(
#'     c(3, 0.05),
#'     c(2.5, 0.05),
#'     c(2, 0.05)
#'   ),
#'   factor_structure = "mixed",
#'   num.factor = "multiple",
#'   real_stats = TRUE,
#'   real_means_vars = list(
#'     c(mean = 5, var = 1),
#'     c(mean = 4.5, var = 0.8),
#'     c(mean = 4.0, var = 0.6)
#'   ),
#'   seed = 123
#' )
#' plot_simData(sim_object = sim2, data = "merged", type = "heatmap")
#'
#' @export
simulateMultiOmics <- function(
    vector_features,
    n_samples,
    n_factors,
    snr = 2,
    signal.samples = c(5, 0.05),
    signal.features = NULL,
    factor_structure = "mixed",
    num.factor = "multiple",
    seed = NULL,
    real_stats = FALSE,
    real_means_vars = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  if (length(vector_features) < 2) stop("Provide at least 2 omics for SUMO simulation")
  k <- length(vector_features)
  if (!is.list(signal.features) || length(signal.features) != k) {
    stop("signal.features must be a list of length equal to vector_features")
  }
  if (!num.factor %in% c("multiple", "single")) stop("num.factor must be either 'multiple' or 'single'")
  if (!factor_structure %in% c("shared", "unique", "mixed", "partial", "custom")) stop("Invalid factor_structure option")
  if (real_stats) {
    if (is.null(real_means_vars) || length(real_means_vars) != k) {
      stop("When real_stats = TRUE, provide real_means_vars = list(c(mean=..., var=...), ...) for each omic.")
    }
  }

  # Sample block assignment
  assigned_indices_samples <- list()
  if (num.factor == "single") {
    block_size <- sample(floor(0.15 * n_samples):ceiling(0.45 * n_samples), 1)
    start_idx <- sample(1:(n_samples - block_size + 1), 1)
    block <- start_idx:(start_idx + block_size - 1)
    assigned_indices_samples <- list(factor1 = block)
    n_factors <- 1
  } else {
    used_indices <- c()
    available_indices <- setdiff(1:n_samples, used_indices)
    i <- 1
    while (i <= n_factors && length(available_indices) > 5) {
      block_size <- sample(5:min(ceiling(n_samples / n_factors) + 5, length(available_indices)), 1)
      possible_starts <- available_indices[available_indices + block_size - 1 <= max(available_indices)]
      if (length(possible_starts) == 0) break
      start_idx <- sample(possible_starts, 1)
      block <- start_idx:(start_idx + block_size - 1)
      if (all(block %in% available_indices)) {
        assigned_indices_samples[[paste0("factor", i)]] <- block
        used_indices <- c(used_indices, block)
        available_indices <- setdiff(1:n_samples, used_indices)
        i <- i + 1
      } else {
        available_indices <- setdiff(available_indices, start_idx)
      }
    }
    while (i <= n_factors && length(available_indices) >= 5) {
      block_size <- sample(5:min(ceiling(n_samples / n_factors) + 5, length(available_indices)), 1)
      possible_starts <- available_indices[available_indices + block_size - 1 <= max(available_indices)]
      if (length(possible_starts) == 0) break
      start_idx <- sample(possible_starts, 1)
      block <- start_idx:(start_idx + block_size - 1)
      if (all(block %in% available_indices)) {
        assigned_indices_samples[[paste0("factor", i)]] <- block
        used_indices <- c(used_indices, block)
        available_indices <- setdiff(1:n_samples, used_indices)
        i <- i + 1
      } else {
        available_indices <- setdiff(available_indices, start_idx)
      }
    }
    while (i <= n_factors) {
      fallback_size <- min(5, length(available_indices))
      if (fallback_size < 5) {
        assigned_indices_samples[[paste0("factor", i)]] <- sample(1:n_samples, 5)
      } else {
        fallback_start <- sample(available_indices[1:(length(available_indices) - fallback_size + 1)], 1)
        fallback_block <- fallback_start:(fallback_start + fallback_size - 1)
        assigned_indices_samples[[paste0("factor", i)]] <- fallback_block
        used_indices <- c(used_indices, fallback_block)
        available_indices <- setdiff(1:n_samples, used_indices)
      }
      i <- i + 1
    }
  }

  # Map factors to omics
  factor_omic_map <- list()
  if (num.factor == "single") {
    if (factor_structure == "shared") {
      factor_omic_map[["factor1"]] <- 1:k
    } else if (factor_structure == "unique") {
      factor_omic_map[["factor1"]] <- sample(1:k, 1)
    } else if (factor_structure == "partial") {
      factor_omic_map[["factor1"]] <- sort(sample(1:k, ifelse(k == 2, 2, sample(2:(k - 1), 1))))
    } else {
      stop("Invalid factor_structure for single factor.")
    }
  } else {
    for (i in 1:n_factors) {
      if (factor_structure == "shared") {
        factor_omic_map[[paste0("factor", i)]] <- 1:k
      } else if (factor_structure == "unique") {
        factor_omic_map[[paste0("factor", i)]] <- sample(1:k, 1)
      } else if (factor_structure == "mixed") {
        factor_omic_map[[paste0("factor", i)]] <- sort(sample(1:k, sample(1:k, 1)))
      } else if (factor_structure == "partial") {
        factor_omic_map[[paste0("factor", i)]] <- sample(1:k, ifelse(k == 2, 2, sample(2:(k - 1), 1)))
      } else {
        stop("Custom factor_structure not supported yet.")
      }
    }
  }

  # Alphas (sample-level latent factors)
  list_alphas <- list()
  for (i in seq_len(n_factors)) {
    alpha_vec <- numeric(n_samples)
    block <- assigned_indices_samples[[paste0("factor", i)]]
    scores_block <- rnorm(length(block), mean = signal.samples[1], sd = signal.samples[2])
    alpha_vec[block] <- scores_block
    non_block <- setdiff(seq_len(n_samples), block)
    alpha_vec[non_block] <- rnorm(length(non_block), mean = 0, sd = 0.05)
    list_alphas[[paste0("alpha", i)]] <- alpha_vec
  }

  # Create omics
  omic.list <- vector("list", k)
  list_betas <- vector("list", k)
  signal_annotation <- list(samples = assigned_indices_samples)

  for (omic_idx in 1:k) {
    n_features <- vector_features[omic_idx]
    feature_mean <- signal.features[[omic_idx]][1]
    feature_sd <- signal.features[[omic_idx]][2]
    omic_data <- matrix(0, nrow = n_samples, ncol = n_features)
    list_betas[[omic_idx]] <- list()
    used_feature_blocks <- list()

    for (factor_i in 1:n_factors) {
      factor_name <- paste0("factor", factor_i)
      if (!(omic_idx %in% factor_omic_map[[factor_name]])) next

      # Feature block helper
      generate_sequential_feature_block <- function(available_indices, min_percent = 0.1, max_percent = 0.15, used_blocks = list()) {
        total_features <- max(available_indices)
        min_block_size <- max(5, ceiling(total_features * min_percent))
        max_block_size <- min(ceiling(total_features * max_percent), length(available_indices))
        block_size <- sample(min_block_size:max_block_size, 1)
        possible_starts <- available_indices[available_indices + block_size - 1 <= max(available_indices)]
        for (start_idx in sample(possible_starts)) {
          block <- start_idx:(start_idx + block_size - 1)
          overlaps <- any(unlist(lapply(used_blocks, function(x) any(block %in% x))))
          if (!overlaps && all(block %in% available_indices)) return(block)
        }
        return(integer(0))
      }

      available_indices <- setdiff(1:n_features, unlist(used_feature_blocks))
      feature_block <- generate_sequential_feature_block(available_indices, used_blocks = used_feature_blocks)
      used_feature_blocks[[factor_name]] <- feature_block
      beta <- rnorm(n_features, mean = 0, sd = 0.01)
      beta[feature_block] <- rnorm(length(feature_block), mean = feature_mean + 0.5 * factor_i, sd = feature_sd)
      list_betas[[omic_idx]][[paste0("beta", factor_i)]] <- beta
      alpha <- list_alphas[[paste0("alpha", factor_i)]]
      omic_data <- omic_data + outer(alpha, beta)
    }

    # Always add noise
    if (real_stats) {
      real_mean <- real_means_vars[[omic_idx]]["mean"]
      real_var <- real_means_vars[[omic_idx]]["var"]
      noise_matrix <- matrix(rnorm(n_samples * n_features, mean = real_mean, sd = sqrt(real_var)), nrow = n_samples, ncol=n_features)
      signal_var <- var(as.vector(omic_data))
      if (signal_var == 0) {
        omic_data <- noise_matrix
      } else {
        scaling_factor <- sqrt((snr * real_var) / signal_var)
        omic_data <- omic_data * scaling_factor + noise_matrix
      }
    } else {
      signal_var <- var(as.vector(omic_data))
      if (signal_var == 0) {
        noise_matrix <- matrix(rnorm(n_samples * n_features, mean = 0, sd = 1), nrow = n_samples)
        omic_data <- noise_matrix
      } else {
        noise_sd <- sqrt(signal_var / snr)
        noise_matrix <- matrix(rnorm(n_samples * n_features, mean = 0, sd = noise_sd), nrow = n_samples)
        omic_data <- omic_data + noise_matrix
      }
    }

    colnames(omic_data) <- paste0("omic", omic_idx, "_feature_", 1:n_features)
    rownames(omic_data) <- paste0("sample_", 1:n_samples)
    omic.list[[omic_idx]] <- omic_data
  }

  signal_annotation$features <- list_betas
  names(omic.list) <- paste0("omic", seq_len(k))
  concatenated_dataset <- do.call(cbind, omic.list)
  rownames(concatenated_dataset) <- paste0("sample_", seq_len(n_samples))

  return(list(
    concatenated_datasets = list(concatenated_dataset),
    omics = omic.list,
    list_alphas = list_alphas,
    list_betas = list_betas,
    signal_annotation = signal_annotation,
    factor_structure = factor_structure,
    factor_map = factor_omic_map
  ))
}
