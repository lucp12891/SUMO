
#' @name feature_selection_one
#' @title Dividing features to create vectors with signal in the first omic
#' @param n_features_one number of features of first omic
#' @param num.factor type of factors - single or multiple
#' @param no_factor number of factors
#'
#'
#' @export
feature_selection_one <- function(n_features_one, num.factor, no_factor) {
  # Ensure the vector has at least 10% elements
  if (n_features_one < 0.1 * n_features_one) {
    stop("The vector must have at least 10% elements.")
  }

  # Create the vector
  vector <- 1:n_features_one

  # Determine the number of elements to select for the first subset
  num_elements <- sample(0.1 * n_features_one:min(0.7 * n_features_one, n_features_one), 1)

  # Select the consecutive elements starting from the first element
  first_subset <- vector[1:num_elements]

  # Handling the division logic for the other subsets
  if (num.factor == 'multiple' && !is.null(no_factor)) {
    num <- no_factor - 1  # Adjust num to account for the first subset that starts with 1
  } else {
    stop("Invalid input for num.factor or no_factor.")
  }

  remaining_length <- n_features_one - length(first_subset)
  if (num * 10 >= remaining_length) {
    stop("Minimum segment size constraint is too large for the given vector length and number of segments.")
  }

  # Generate breakpoints for the remaining length
  breakpoints <- if (num > 1) {
    sort(sample(1:(remaining_length - num * 10), num - 1, replace = FALSE)) + (0:(num - 2)) * 10
  } else {
    numeric(0)
  }

  # Calculate the sizes of each segment
  segment_sizes <- diff(c(0, breakpoints, remaining_length))

  # Ensure all segment sizes are greater than the minimum size
  if (any(segment_sizes <= 10)) {
    return(feature_selection_one(n_features_one, num.factor, no_factor))
  }

  # Calculate start and end points of segments
  end_points <- cumsum(segment_sizes) + length(first_subset)
  start_points <- c(length(first_subset) + 1, end_points[-length(end_points)] + 1)

  # Generate the full vector starting after the first subset
  full_vector <- vector[(length(first_subset) + 1):n_features_one]

  # Extract sub-vectors based on start and end points
  vectors <- mapply(function(start, end) {
    full_vector[(start - length(first_subset)):(end - length(first_subset))]
  }, start_points, end_points, SIMPLIFY = FALSE)

  # Create a function to select 40% from each segment
  select_40_percent_vector <- function(original_vector) {
    target_length <- ceiling(0.4 * length(original_vector))  # Selecting 40% of the original vector length
    start_index <- sample(1:(length(original_vector) - target_length + 1), 1)  # Randomly select the start index
    original_vector[start_index:(start_index + target_length - 1)]  # Extract the sub-vector
  }

  # Create the final segments
  sub_vectors <- list(first_subset)

  for (segment in vectors) {
    sub_vectors <- append(sub_vectors, list(select_40_percent_vector(segment)))
  }

  return(sub_vectors)
}

#no_factor = 4
#tet = feature_selection_one(n_features_one = 5000, num.factor = "multiple", no_factor)
# Check if there is an intersection between any pair of vectors in the list
has_intersection <- function(lst) {
  # Compare each pair of vectors to see if there is an intersection
  for (i in 1:(length(lst) - 1)) {
    for (j in (i + 1):length(lst)) {
      # Check intersection between vec i and vec j
      if (length(intersect(lst[[i]], lst[[j]])) > 0) {
        return(TRUE) # Return TRUE if there is an intersection
      }
    }
  }
  return(FALSE) # Return FALSE if no intersections found
}
# Check for intersections in the list
#result <- has_intersection(tet)

# Print the result
#print(result)


#feature_selection_one(n_features_one = 5000, num.factor = "multiple", no_factor  = 3)
