#' @name simulateMultiOmics
#' @title Simulation of omics with predefined single or multiple latent factors in multi-omics
#'
#' @description Simulate multiple omics (>2) datasets with non-overlapping sample and feature signal blocks.
#'
#' @details This function generates synthetic omics data where each omic layer has its own feature space and noise characteristics. The sample signal blocks for each latent factor are non-overlapping and sequential with random gaps. Feature signal blocks are generated per omic with sequential non-overlapping segments.
#'
#' @param vector_features Integer vector indicating number of features per omic (length k for k omics).
#' @param n_samples Total number of samples across all omics.
#' @param n_factors Number of latent factors.
#' @param snr Signal-to-noise ratio.
#' @param signal.samples Mean and SD for generating sample signal values (e.g., c(mean, sd)).
#' @param signal.features List of vectors with mean and SD for features per omic (e.g., list(c(3,0.2), c(2.5,0.15))).
#' @param factor_structure Character. "shared", "exclusive", "mixed", "partial", or "custom" factor distribution
#' @param num.factor Character. "multiple" (default) or "single"
#' @param seed Optional. Set random seed for reproducibility.
#' @importFrom utils head
#' @importFrom stats var
#' @importFrom stats rnorm setNames runif
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom gridExtra grid.arrange
#' @importFrom rlang abort
#' @importFrom stats sd
#' @importFrom utils data
#' @importFrom readxl read_excel
#' @importFrom readr read_csv
#' @importFrom stringr str_detect
#' @include divide_samples.R
#' @include feature_selection_one.R
#' @include feature_selection_two.R
#' @include divide_features_one.R
#' @include divide_features_two.R
#' @return A list containing:
#' \itemize{
#'   \item \code{omic.list}: List of simulated omic datasets.
#'   \item \code{signal_annotation}: Signal sample indices per factor.
#'   \item \code{list_alphas}, \code{list_betas}: Latent factor loading vectors.
#' }
#'
#' @examples
#' sim_object1 <- simulateMultiOmics(
#'   vector_features = c(3000, 2500, 2000),
#'   n_samples = 100,
#'   n_factors = 3,
#'   snr = 3,
#'   signal.samples = c(5, 1),
#'   signal.features = list(
#'     c(3, 0.3),   # omic1 signal mean/sd
#'     c(2.5, 0.25),# omic2 signal mean/sd
#'     c(2, 0.2)    # omic3 signal mean/sd
#'   ),
#'   factor_structure = "mixed",
#'   num.factor = "multiple",
#'   seed = 123
#' )
#'
#' # View available elements
#' names(sim_object1)
#'
#' # Visualize the simulated data
#' plot_simData(sim_object = sim_object1, data = "merged", type = "heatmap")
#'
#' sim_object2 <- simulateMultiOmics(
#'   vector_features = c(3000, 2500),
#'   n_samples = 100,
#'   n_factors = 1,
#'   snr = 0.5,
#'   signal.samples = c(3, 1),
#'   signal.features = list(
#'     c(3.5, 0.3),   # omic1 signal mean/sd
#'     c(4, 0.2)    # omic3 signal mean/sd
#'   ),
#'   factor_structure = "shared",
#'   num.factor = "single",
#'   seed = NULL
#' )
#'
#' # Visualize the simulated data
#' plot_simData(sim_object = sim_object2, data = "merged", type = "heatmap")
#'
#' @export
simulateMultiOmics <- function(
    vector_features,
    n_samples,
    n_factors,
    snr = 2,
    signal.samples = c(5, 1),
    signal.features = NULL,
    factor_structure = "mixed",
    num.factor = "multiple",
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Validate input
  if (length(vector_features) < 2) stop("Provide at least 2 omics for SUMO simulation")
  k <- length(vector_features)
  if (!is.list(signal.features) || length(signal.features) != k) {
    stop("signal.features must be a list of length equal to vector_features")
  }
  if (!num.factor %in% c("multiple", "single")) stop("num.factor must be either 'multiple' or 'single'")
  if (!factor_structure %in% c("shared", "unique", "mixed", "partial", "custom")) stop("Invalid factor_structure option")

  # Step 1: Divide samples into non-overlapping blocks for each factor with sequential indices and random gaps
  if (num.factor == "single") {
    assigned_indices_samples <- list()
    block_size <- sample(
      floor(0.15 * n_samples):ceiling(0.45 * n_samples),
      1
    )
    start_idx <- sample(1:(n_samples - block_size + 1), 1)
    block <- start_idx:(start_idx + block_size - 1)

    assigned_indices_samples <- list(factor1 = block)
    n_factors <- 1  # Override to 1 factor
  } else {
    assigned_indices_samples <- list()
    used_indices <- c()
    max_range <- n_samples
    available_indices <- setdiff(1:max_range, used_indices)

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
        available_indices <- setdiff(1:max_range, used_indices)
        i <- i + 1
      } else {
        # skip block if overlap found
        available_indices <- setdiff(available_indices, start_idx)
      }
    }

    # Attempt to assign empty factor slots with remaining space if needed
    while (i <= n_factors && length(available_indices) >= 5) {
      block_size <- sample(5:min(ceiling(n_samples / n_factors) + 5, length(available_indices)), 1)
      possible_starts <- available_indices[available_indices + block_size - 1 <= max(available_indices)]

      if (length(possible_starts) == 0) break

      start_idx <- sample(possible_starts, 1)
      block <- start_idx:(start_idx + block_size - 1)

      if (all(block %in% available_indices)) {
        assigned_indices_samples[[paste0("factor", i)]] <- block
        used_indices <- c(used_indices, block)
        available_indices <- setdiff(1:max_range, used_indices)
        i <- i + 1
      } else {
        available_indices <- setdiff(available_indices, start_idx)
      }
    }

    # Final fallback to fill any remaining factor slots with random valid blocks
    while (i <= n_factors) {
      fallback_size <- min(5, length(available_indices))
      if (fallback_size < 5) {
        assigned_indices_samples[[paste0("factor", i)]] <- sample(1:n_samples, 5)
      } else {
        fallback_start <- sample(available_indices[1:(length(available_indices) - fallback_size + 1)], 1)
        fallback_block <- fallback_start:(fallback_start + fallback_size - 1)
        assigned_indices_samples[[paste0("factor", i)]] <- fallback_block
        used_indices <- c(used_indices, fallback_block)
        available_indices <- setdiff(1:max_range, used_indices)
      }
      i <- i + 1
    }
  }
  # Initialize tracking containers
  omic.list <- vector("list", k)
  list_alphas <- list()
  list_betas <- vector("list", k)
  signal_annotation <- list(samples = assigned_indices_samples)

  # Continue with latent factor design...
  # (to be implemented in next chunk)

  # Step 2: Define factor allocation across omics
  factor_omic_map <- list()

  if (num.factor == "single") {
    if (factor_structure == "shared") {
      factor_omic_map[["factor1"]] <- 1:k
    } else if (factor_structure == "unique") {
      factor_omic_map[["factor1"]] <- sample(1:k, 1)
    } else if (factor_structure == "partial") {
      if(k == 2){
        n_partial = 2
        factor_omic_map[["factor1"]] <- sort(sample(1:k, n_partial))
      }else{
        n_partial <- sample(2:(k - 1), 1)
        factor_omic_map[["factor1"]] <- sort(sample(1:k, n_partial))
      }
    } else {
      stop("Invalid factor_structure for single factor. Choose from 'shared', 'unique', or 'partial'.")
    }
  } else {
    for (i in 1:n_factors) {
      if (factor_structure == "shared") {
        factor_omic_map[[paste0("factor", i)]] <- 1:k
      } else if (factor_structure == "unique") {
        assigned <- sample(1:k, 1)
        factor_omic_map[[paste0("factor", i)]] <- assigned
      } else if (factor_structure == "mixed") {
        assigned <- sort(sample(1:k, sample(1:k, 1)))
        factor_omic_map[[paste0("factor", i)]] <- assigned
      } else if (factor_structure == "partial") {
        if(k == 2){
          assigned = 2
          factor_omic_map[[paste0("factor", i)]] <- assigned
        }else{
          assigned <- sample(2:(k - 1), 1)
          factor_omic_map[[paste0("factor", i)]] <- assigned
        }
      } else if (factor_structure == "custom") {
        stop("Custom factor_structure currently not supported in this version")
      }
    }
  }

  # Step 3: Create alphas (sample scores) for each factor
  list_alphas <- list()
  if (num.factor == "single") {
    alpha_vec <- numeric(n_samples)
    block <- assigned_indices_samples[["factor1"]]
    scores <- rnorm(length(block), mean = signal.samples[1], sd = signal.samples[2])
    alpha_vec[block] <- scores
    list_alphas[["alpha1"]] <- alpha_vec
  } else {
    for (i in seq_len(n_factors)) {
      alpha_vec <- numeric(n_samples)
      block <- assigned_indices_samples[[paste0("factor", i)]]

      # Assign signal scores to block members
      scores_block <- rnorm(length(block), mean = signal.samples[1], sd = signal.samples[2])
      alpha_vec[block] <- scores_block

      # Assign small-noise scores to non-block members
      non_block <- setdiff(seq_len(n_samples), block)
      alpha_vec[non_block] <- rnorm(length(non_block), mean = 0, sd = 0.05)

      list_alphas[[paste0("alpha", i)]] <- alpha_vec
    }
  }

  # Step 4: Generate feature weights (betas) and construct omic-specific datasets
  for (omic_idx in 1:k) {
    n_features <- vector_features[omic_idx]
    feature_mean <- signal.features[[omic_idx]][1]
    feature_sd <- signal.features[[omic_idx]][2]

    omic_data <- matrix(0, nrow = n_samples, ncol = n_features)
    list_betas[[omic_idx]] <- list()

    # Keep track of used feature blocks for this omic
    used_feature_blocks <- list()

    for (factor_i in 1:n_factors) {

      # sequential function
      generate_sequential_feature_block <- function(available_indices, min_percent = 0.1, max_percent = 0.15, used_blocks = list()) {
        total_features <- max(available_indices)
        min_block_size <- max(5, ceiling(total_features * min_percent))
        max_block_size <- min(ceiling(total_features * max_percent), length(available_indices))

        block_size <- sample(min_block_size:max_block_size, 1)
        possible_starts <- available_indices[available_indices + block_size - 1 <= max(available_indices)]

        for (start_idx in sample(possible_starts)) {
          block <- start_idx:(start_idx + block_size - 1)

          # Check overlap
          overlaps <- any(unlist(lapply(used_blocks, function(x) any(block %in% x))))
          if (!overlaps && all(block %in% available_indices)) {
            return(block)
          }
        }

        return(integer(0))  # Return empty if no valid block
      }
      # end function

      factor_name <- paste0("factor", factor_i)
      if (!(omic_idx %in% factor_omic_map[[factor_name]])) next

      # Generate a non-overlapping, sequential feature block
      available_indices <- setdiff(1:n_features, unlist(used_feature_blocks))
      feature_block <- generate_sequential_feature_block(available_indices, used_blocks = used_feature_blocks)

      # Store for future overlap checks
      used_feature_blocks[[factor_name]] <- feature_block

      # Assign feature weights
      beta <- rnorm(n_features, mean = 0, sd = 0.01)
      beta[feature_block] <- rnorm(length(feature_block), mean = feature_mean + 0.5 * factor_i, sd = feature_sd)
      list_betas[[omic_idx]][[paste0("beta", factor_i)]] <- beta

      alpha <- list_alphas[[paste0("alpha", factor_i)]]
      omic_data <- omic_data + outer(alpha, beta)
    }

    # Add noise
    signal_variance <- var(as.vector(omic_data))
    noise_sd <- sqrt(signal_variance / snr)
    noise_matrix <- matrix(rnorm(n_samples * n_features, mean = 0, sd = noise_sd), nrow = n_samples)
    omic_data <- omic_data + noise_matrix

    # Store
    colnames(omic_data) <- paste0("omic", omic_idx, "_feature_", 1:n_features)
    rownames(omic_data) <- paste0("sample_", 1:n_samples)
    omic.list[[omic_idx]] <- omic_data
  }

  signal_annotation$features <- list_betas

  # Step 5: Finalize simulated object
  names(omic.list) <- paste0("omic", seq_len(k))

  concatenated_dataset <- do.call(cbind, omic.list)
  rownames(concatenated_dataset) <- paste0("sample_", seq_len(n_samples))

  sim_object <- list(
    concatenated_datasets = list(concatenated_dataset),
    omics = omic.list,
    list_alphas = list_alphas,
    list_betas = list_betas,
    signal_annotation = signal_annotation,
    factor_structure = factor_structure,
    factor_map = factor_omic_map
  )

  return(sim_object)

}

# Load the M_SUMO package (assumes it's loaded or sourced already)

# Simulate 3 omics datasets with:
# - 3000, 2500, and 2000 features respectively
# - 100 samples
# - 4 latent factors
# - Mixed factor structure
# - Moderate signal-to-noise ratio (SNR = 3)
