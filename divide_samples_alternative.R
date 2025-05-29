#' Divide Sample Indices into Subsets with Minimum Size Constraints
#'
#' This utility function divides a sequence of sample indices into `num` segments
#' ensuring that each segment meets a specified minimum size. It optionally
#' extracts a subset of each segment based on predefined selection logic:
#' - For a single group (`num = 1`): selects a random contiguous sub-vector
#'   comprising between 10% and 55% of the total samples.
#' - For multiple groups (`num > 1`): selects a contiguous sub-vector
#'   comprising approximately 75% of each segment.
#'
#' This function is primarily used for randomized simulation of sample blocks,
#' useful in bootstrapping, subsampling, or simulating latent factor scores
#' across multi-omics datasets.
#'
#' @param n_samples Integer. Total number of samples to divide.
#' @param num Integer. Number of desired segments or latent factors.
#' @param min_size Integer. Minimum size (length) allowed for each segment.
#'
#' @return A list of integer vectors. Each vector contains a sequence of indices
#' representing a subsample of the corresponding segment.
#'
#' @examples
#' divide_samples(n_samples = 100, num = 3, min_size = 10)
#' divide_samples(n_samples = 50, num = 1, min_size = 5)
#'
#' @export
divide_samples <- function(n_samples, num, min_size) {
  # Validate input
  if (!is.numeric(n_samples) || !is.numeric(num) || !is.numeric(min_size)) {
    stop("All inputs must be numeric.")
  }
  if (num * min_size >= n_samples) {
    stop("Minimum segment size constraint is too large for the given number of samples and segments.")
  }

  # Helper for single vector scenario
  if (num == 1) {
    vector <- 1:n_samples
    min_required <- ceiling(0.1 * n_samples)
    if (length(vector) < min_required) {
      stop("The vector must have at least 10% of the total elements.")
    }
    # Sample size between 10% and 55% of the full vector
    num_elements <- sample(min_required:min(ceiling(0.55 * n_samples), length(vector)), 1)
    start_index <- sample(1:(length(vector) - num_elements + 1), 1)
    sub_vector <- vector[start_index:(start_index + num_elements - 1)]
    return(sub_vector)
  }

  # Proceed for num > 1
  remaining_length <- n_samples - num * min_size
  breakpoints <- if (remaining_length > 0) {
    sort(sample(1:remaining_length, num - 1, replace = FALSE))
  } else {
    numeric(num - 1)
  }

  # Adjust breakpoints by adding min_size offsets
  breakpoints <- breakpoints + (0:(num - 2)) * min_size

  # Determine segment sizes
  segment_sizes <- diff(c(0, breakpoints, n_samples))

  if (any(segment_sizes <= min_size)) {
    return(divide_samples(n_samples, num, min_size)) # Retry with recursion
  }

  end_points <- cumsum(segment_sizes)
  start_points <- c(1, head(end_points, -1) + 1)
  full_vector <- 1:n_samples

  # Generate segments
  segments <- Map(function(start, end) full_vector[start:end], start_points, end_points)

  # Extract 75% sub-vectors from each segment
  select_subvector <- function(vec, proportion = 0.75) {
    target_length <- ceiling(proportion * length(vec))
    start_index <- sample(1:(length(vec) - target_length + 1), 1)
    vec[start_index:(start_index + target_length - 1)]
  }

  sub_vectors <- lapply(segments, select_subvector)
  return(sub_vectors)
}

#divide_samples(n_samples = 100, num = 8, min_size = 10)
