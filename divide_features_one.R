#' @name divide_features_one
#' @title Dividing features to create vectors with signal in the first omic for single data
#' @param n_features_one number of features of first omic
#' @param num.factor number of factor = '1'
#' @include divide_vector.R
#'
#' @export
divide_features_one <- function(n_features_one, num.factor) {
  min_size = 10
  if (num.factor == 'single'){
    num = 1
  }else if(num.factor == 'multiple'){
    num = n_factors
  }
  if (num * min_size >= n_features_one) {
    stop("Minimum segment size constraint is too large for the given vector length and number of segments.")
  }

  # Initialize the breakpoints
  breakpoints <- numeric(num - 1)

  # Calculate the remaining length to be divided after accounting for minimum sizes
  remaining_length <- n_features_one - num * min_size

  # Generate breakpoints for the remaining length
  if (remaining_length > 0) {
    breakpoints <- sort(sample(1:remaining_length, num - 1, replace = FALSE))
  }

  # Adjust breakpoints to account for minimum sizes
  breakpoints <- breakpoints + (0:(num - 2)) * min_size

  # Calculate the sizes of each segment
  segment_sizes <- diff(c(0, breakpoints, n_features_one))

  # Ensure all segment sizes are greater than min_size
  if (any(segment_sizes <= min_size)) {
    return(divide_vector(n_features_one, num, min_size)) # Retry if any segment is less than or equal to min_size
  }

  # Calculate the cumulative sum to get end points of segments
  end_points <- cumsum(segment_sizes)
  start_points <- c(1, end_points[-length(end_points)] + 1)

  # Generate the full vector
  full_vector <- 1:n_features_one

  # Extract sub-vectors based on start and end points
  vectors <- mapply(function(start, end) {
    full_vector[start:end]
  }, start_points, end_points, SIMPLIFY = FALSE)

  #return(sub_vectors)
  sub_vectors <- list()

  select_80_percent_vector <- function(original_vector) {
    target_length <- ceiling(0.8 * length(original_vector)) # Selecting 80% of the original vector length
    start_index <- sample(1:(length(original_vector)-target_length), 1)  # Randomly select the start index
    # select the 80% of the original vector
    sub_vector <- original_vector[start_index:(start_index + target_length - 1)]  # Extract the sub-vector
    return(sub_vector)
  }

  # creating the segments
  for(i in seq_along(vectors)){

    segment <- vectors[[i]]
    sub_vector <- select_80_percent_vector(segment)
    sub_vectors[[i]] <- sub_vector
  }

  return(sub_vectors)
}
# Suppressing the global variable warning
utils::globalVariables(c("n_factors"))
