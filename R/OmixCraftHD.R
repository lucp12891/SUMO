#' @name OmixCraftHD
#' @title Simulation of high-dimensional data with predefined single factor or multiple factors in multi-omics
#' @param vector_features Vector of features assigned to the two simulated datasets respectively '1' first dataset, '2' second dataset
#' @param n_samples The number of samples common between the two simulated datasets
#' @param sigmas_vector Vector for the noise variability for the two simulated datasets respectively,  '1' first dataset, '2' second dataset
#' @param n_factors Number of predefined factors
#' @param num.factor Category of factors to be simulated specified as 'single', or 'multiple'.
#' @param advanced_dist Applicable only when num.factor = 'multiple'. Contains six possible arguments, '', NULL, 'mixed', 'omic.one', or 'omic.two', 'exclusive'
#' @importFrom stats rnorm setNames runif
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom gridExtra grid.arrange
#' @importFrom rlang abort
#' @include divide_samples.R
#' @include feature_selection_one.R
#' @include feature_selection_two.R
#' @include divide_features_one.R
#' @include divide_features_two.R
#' @examples
#' # Examples
#' set.seed(1234)
#' output_obj <- OmixCraftHD(
#'   vector_features = c(2000,3000),
#'   sigmas_vector=c(8,5),
#'   n_samples=100,
#'   n_factors=5,
#'   num.factor='multiple',
#'   advanced_dist='mixed'
#' )
#' output_obj <- OmixCraftHD(
#'   vector_features = c(5000,3000),
#'   sigmas_vector=c(3,4),
#'   n_samples=30, n_factors=1
#' )
#' @export
# @param type.factor Applicable only when num.factor = 'single'. type of factor needed to be simulated. Contains two type 'shared', 'unique'. 'shared' refers to latent factor present in both the dataset. 'unique' refers to latent factor present in one of the datasets.
# @param signal_location Applicable only when num.factor = 'single' when type.factor = 'unique'. Contains three possible arguments empty (''), omic.one'  or 'omic.two' with 'omic.one' refers to factor only in the first data while 'omic.two' indicates the factor present only on the second data.

OmixCraftHD <- function(vector_features = c(2000,2000), n_samples = 50, sigmas_vector = c(3,5), n_factors = 3, num.factor = 'multiple', advanced_dist = NULL){
  if(num.factor == 'single'){
    # FeaturesD
    n_features_one <- vector_features[1]
    n_features_two <- vector_features[2]

    n_samples = n_samples

    # Samples scores
    num = 1
    #min_size = ceiling(n_samples*0.1)
    #divide_samples <- divide_samples(n_samples, num, min_size)
    select_consecutive_samples <- function(n_samples) {
      # Ensure the vector has at least 10% elements
      vector = 1:n_samples
      if (length(vector) < ceiling(0.1*n_samples)) {
        stop("The vector must have at least 10% elements.")
      }

      # Determine the number of elements to select (between 10% and 70% of the elements)
      num_elements <- sample(ceiling(0.1*n_samples):min(ceiling(0.7*n_samples), length(vector)), 1)

      # Determine the starting point for the consecutive elements
      start_index <- sample(1:(length(vector) - num_elements + 1), 1)

      # Select the consecutive elements
      selected_elements <- vector[start_index:(start_index + num_elements - 1)]

      return(selected_elements)
    }

    # Example usage:
    divide_samples <- select_consecutive_samples(n_samples=n_samples)
    assigned_indices_samples = divide_samples
    #Explore features
    n_factors = 1

    type.factor <- readline(prompt = "Please provide type factor; unique or shared or NULL: ")

    if (type.factor == 'unique'){

      signal_location <- readline(prompt = "Please provide signal location; omic.one or omic.two: ")
      #sigma <- readline(prompt = "Please provide variance for the backgound noise (number): ")

      if(signal_location == ""){

        omic <- c('omic.one', 'omic.two')
        selected.omic <- sample(omic, 1)

        if (selected.omic == 'omic.one'){
          alpha <- list()
          for (i in seq_along(assigned_indices_samples)) {
            alpha[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
          }

          # Assigning the signal to the selected indices
          for (i in seq_along(alpha)) {
            indices <- assigned_indices_samples#[[i]]
            alpha[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
          }

          # # First OMIC data
          # FeaturesD signal start_and_end
          assigned_indices_features = divide_features_one(n_features_one, num.factor = 1)    # feature indices

          all_indices <- seq_along(assigned_indices_features)

          # Create only one vector in list_betas using the random index
          beta <- list()
          for (i in all_indices) {
            beta[[paste0("beta", i)]] <- rnorm(n_features_one, 0, 0.05)
          }

          # Assign corresponding values to betas variables based on assigned_indices_features
          for (i in all_indices) {
            indices <- assigned_indices_features[[i]]
            beta[[i]][indices] <- rnorm(length(indices), (4.0 + 0.5 * i), 0.05)  # Adjust values dynamically
          }

          pattern_alpha <-"^alpha\\d+"
          matches_alpha <- grepl(pattern_alpha, names(alpha), perl = TRUE)
          n_alpha = length(matches_alpha[matches_alpha == TRUE])

          pattern_beta <- "^beta\\d+"
          matches_beta <- grepl(pattern_beta, names(beta), perl = TRUE)
          n_beta = length(matches_beta[matches_beta == TRUE])

          # Initialize the data list to store the results
          data_list_i <- list()

          # Loop through each alpha and beta combination
          for (i in 1:min(length(alpha), length(beta))) {
            data_i <- alpha[[i]] %*% t(beta[[i]][1:n_features_one])
            data_list_i[[paste0("data.", i)]] <- data_i
          }

          # Combine the results into a single data variable
          data.1 <- Reduce(`+`, data_list_i)

          omic.one <- list()

          # Noise in the first detasets
          sigma <- sigmas_vector[1]

          eps1 <- rnorm(n_samples*n_features_one, 0, sigma) # noise
          omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
          colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
          rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

          # Second OMIC
          # Generate random gamma values based on the max_factors
          gamma <- list()
          gamma <- alpha #list(), assign the values of alpha for the second data

          assigned_indices_features_omic.two <- divide_features_two(n_features_two, num.factor = 1)

          all_indices <- seq_along(assigned_indices_features_omic.two)

          # Create only one vector in list_betas using the random index
          delta <- list()
          for (i in all_indices) {
            delta[[paste0("delta", i)]] <- rnorm(n_features_two, 0, 0.05)
          }

          # Assign corresponding values to betas variables based on assigned_indices_features
          for (i in all_indices) {
            indices <- assigned_indices_features_omic.two[[i]]
            delta[[i]][indices] <- rnorm(length(indices), 0, 0.05)  # Adjust values dynamically
          }

          pattern_gamma <-"^alpha\\d+"
          matches_gamma <- grepl(pattern_gamma, names(gamma), perl = TRUE)
          n_gamma = length(matches_gamma[matches_gamma == TRUE])

          pattern_delta <- "^delta\\d+"
          matches_delta <- grepl(pattern_delta, names(delta), perl = TRUE)
          n_delta = length(matches_delta[matches_delta == TRUE])

          # Initialize the data list to store the results
          data_list_j <- list()

          # Loop through each alpha and beta combination
          for (j in 1:min(length(gamma), length(delta))) {
            data_j <- gamma[[j]] %*% t(delta[[j]][1:n_features_two])
            data_list_j[[paste0("data.", j)]] <- data_j
          }

          # Combine the results into a single data variable
          data.2 <- Reduce(`+`, data_list_j)

          omic.two <- list()

          # Noise in the second detasets
          sigma <- sigmas_vector[2]

          eps2 <- rnorm(n_samples*n_features_two, 0, sigma) # noise
          omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
          colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
          rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

          phi <- ceiling(n_features_one)

          simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

          concatenated_datasets <- list()
          for (i in 1:length(simulated_datasets$object.two)) {
            concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
            concatenated_datasets[[i]] <- concatenated_data
          }
          sim_output <- return(list(concatenated_datasets = concatenated_datasets, divide_samples = divide_samples, assigned_indices_features=assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, list_alphas=alpha, list_gammas=gamma, list_betas=beta, list_deltas=delta))
          return(sim_output)

        }else if (selected.omic == 'omic.two') {
          alpha <- list()
          for (i in seq_along(assigned_indices_samples)) {
            alpha[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
          }

          # Assigning the signal to the selected indices
          for (i in seq_along(alpha)) {
            indices <- assigned_indices_samples#[[i]]
            alpha[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
          }


          assigned_indices_features = divide_features_one(n_features_one, num.factor = 1)

          all_indices <- seq_along(assigned_indices_features)

          # Create only one vector in list_betas using the random index
          beta <- list()
          for (i in all_indices) {
            beta[[paste0("beta", i)]] <- rnorm(n_features_one, 0, 0.05)
          }

          # Assign corresponding values to betas variables based on assigned_indices_features
          for (i in all_indices) {
            indices <- assigned_indices_features[[i]]
            beta[[i]][indices] <- rnorm(length(indices), 0, 0.05)  # Adjust values dynamically
          }

          pattern_alpha <-"^alpha\\d+"
          matches_alpha <- grepl(pattern_alpha, names(alpha), perl = TRUE)
          n_alpha = length(matches_alpha[matches_alpha == TRUE])

          pattern_beta <- "^beta\\d+"
          matches_beta <- grepl(pattern_beta, names(beta), perl = TRUE)
          n_beta = length(matches_beta[matches_beta == TRUE])

          # Initialize the data list to store the results
          data_list_i <- list()

          # Loop through each alpha and beta combination
          for (i in 1:min(length(alpha), length(beta))) {
            data_i <- alpha[[i]] %*% t(beta[[i]][1:n_features_one])
            data_list_i[[paste0("data.", i)]] <- data_i
          }

          # Combine the results into a single data variable
          data.1 <- Reduce(`+`, data_list_i)

          omic.one <- list()

          # Noise in the first detasets
          sigma <- sigmas_vector[1]

          eps1 <- rnorm(n_samples*n_features_one, 0, sigma) # noise
          omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
          colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
          rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

          gamma <- list()
          gamma <- alpha #list(), assign the values of alpha for the second data

          assigned_indices_features_omic.two <- divide_features_two(n_features_two, num.factor = 1)

          #all_indices <- seq_along(assigned_indices_features_omic.two)

          # Create only one vector in list_betas using the random index
          delta <- list()
          for (i in seq_along(num)) {
            delta[[paste0("delta", i)]] <- rnorm(n_features_two, 0, 0.05)
          }

          # Assign corresponding values to betas variables based on assigned_indices_features
          for (i in seq_along(num)) {
            indices <- assigned_indices_features_omic.two[[i]]
            delta[[i]][indices] <- rnorm(length(indices), 9 + 0.05*i, 0.05)  # Adjust values dynamically
          }

          pattern_gamma <-"^alpha\\d+"
          matches_gamma <- grepl(pattern_gamma, names(gamma), perl = TRUE)
          n_gamma = length(matches_gamma[matches_gamma == TRUE])

          pattern_delta <- "^delta\\d+"
          matches_delta <- grepl(pattern_delta, names(delta), perl = TRUE)
          n_delta = length(matches_delta[matches_delta == TRUE])

          # Initialize the data list to store the results
          data_list_j <- list()

          # Loop through each alpha and beta combination
          for (j in 1:min(length(gamma), length(delta))) {
            data_j <- gamma[[j]] %*% t(delta[[j]][1:n_features_two])
            data_list_j[[paste0("data.", j)]] <- data_j
          }

          # Combine the results into a single data variable
          data.2 <- Reduce(`+`, data_list_j)

          omic.two <- list()

          # Noise in the second detasets
          sigma <- sigmas_vector[2]

          eps2 <- rnorm(n_samples*n_features_two, 0, sigma) # noise
          omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
          colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
          rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

          simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

          concatenated_datasets <- list()
          for (i in 1:length(simulated_datasets$object.two)) {
            concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
            concatenated_datasets[[i]] <- concatenated_data
          }
          sim_output <- return(list(concatenated_datasets = concatenated_datasets, divide_samples = divide_samples, assigned_indices_features=assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, list_alphas=alpha, list_gammas=gamma, list_betas=beta, list_deltas=delta))
          return(sim_output)

        } else{
          stop("Location must be either 'omic.one', 'omic.two' or 'NULL'.")
        }

      }else if (signal_location == 'omic.one'){
        alpha <- list()
        for (i in seq_along(num)) {
          alpha[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
        }

        # Assigning the signal to the selected indices
        for (i in seq_along(alpha)) {
          indices <- assigned_indices_samples
          alpha[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
        }

        # # First OMIC data
        # FeaturesD signal start_and_end
        assigned_indices_features = divide_features_one(n_features_one, num.factor = 1)    # feature indices

        all_indices <- seq_along(assigned_indices_features)

        # Create only one vector in list_betas using the random index
        beta <- list()
        for (i in all_indices) {
          beta[[paste0("beta", i)]] <- rnorm(n_features_one, 0, 0.05)
        }

        # Assign corresponding values to betas variables based on assigned_indices_features
        for (i in all_indices) {
          indices <- assigned_indices_features[[i]]
          beta[[i]][indices] <- rnorm(length(indices), (4.0 + 0.5 * i), 0.05)  # Adjust values dynamically
        }

        pattern_alpha <-"^alpha\\d+"
        matches_alpha <- grepl(pattern_alpha, names(alpha), perl = TRUE)
        n_alpha = length(matches_alpha[matches_alpha == TRUE])

        pattern_beta <- "^beta\\d+"
        matches_beta <- grepl(pattern_beta, names(beta), perl = TRUE)
        n_beta = length(matches_beta[matches_beta == TRUE])

        # Initialize the data list to store the results
        data_list_i <- list()

        # Loop through each alpha and beta combination
        for (i in 1:min(length(alpha), length(beta))) {
          data_i <- alpha[[i]] %*% t(beta[[i]][1:n_features_one])
          data_list_i[[paste0("data.", i)]] <- data_i
        }

        # Combine the results into a single data variable
        data.1 <- Reduce(`+`, data_list_i)

        omic.one <- list()

        # Noise in the first detasets
        sigma <- sigmas_vector[1]

        eps1 <- rnorm(n_samples*n_features_one, 0, sigma) # noise
        omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
        colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
        rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

        # Second OMIC
        # Generate random gamma values based on the max_factors
        gamma <- list()
        gamma <- alpha #list(), assign the values of alpha for the second data

        assigned_indices_features_omic.two <- divide_features_two(n_features_two, num.factor = 1)

        all_indices <- seq_along(assigned_indices_features_omic.two)

        # Create only one vector in list_betas using the random index
        delta <- list()
        for (i in seq_along(num)) {
          delta[[paste0("delta", i)]] <- rnorm(n_features_two, 0, 0.05)
        }

        # Assign corresponding values to betas variables based on assigned_indices_features
        for (i in seq_along(num)) {
          indices <- assigned_indices_features_omic.two[[i]]
          delta[[i]][indices] <- rnorm(length(indices), 0, 0.05)  # Adjust values dynamically
        }

        pattern_gamma <-"^alpha\\d+"
        matches_gamma <- grepl(pattern_gamma, names(gamma), perl = TRUE)
        n_gamma = length(matches_gamma[matches_gamma == TRUE])

        pattern_delta <- "^delta\\d+"
        matches_delta <- grepl(pattern_delta, names(delta), perl = TRUE)
        n_delta = length(matches_delta[matches_delta == TRUE])

        # Initialize the data list to store the results
        data_list_j <- list()

        # Loop through each alpha and beta combination
        for (j in 1:min(length(gamma), length(delta))) {
          data_j <- gamma[[j]] %*% t(delta[[j]][1:n_features_two])
          data_list_j[[paste0("data.", j)]] <- data_j
        }

        # Combine the results into a single data variable
        data.2 <- Reduce(`+`, data_list_j)

        omic.two <- list()

        # Noise in the second detasets
        sigma <- sigmas_vector[2]

        eps2 <- rnorm(n_samples*n_features_two, 0, sigma) # noise
        omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
        colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
        rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

        phi <- ceiling(n_features_one)

        simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

        concatenated_datasets <- list()
        for (i in 1:length(simulated_datasets$object.two)) {
          concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
          concatenated_datasets[[i]] <- concatenated_data
        }
        sim_output <- return(list(concatenated_datasets = concatenated_datasets, divide_samples = divide_samples, assigned_indices_features=assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, list_alphas=alpha, list_gammas=gamma, list_betas=beta, list_deltas=delta))
        return(sim_output)

      }else if (signal_location == 'omic.two') {

        alpha <- list()
        for (i in seq_along(num)) {
          alpha[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
        }

        # Assigning the signal to the selected indices
        for (i in seq_along(alpha)) {
          indices <- assigned_indices_samples#[[i]]
          alpha[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
        }


        assigned_indices_features = divide_features_one(n_features_one, num.factor = 1)

        all_indices <- seq_along(assigned_indices_features)

        # Create only one vector in list_betas using the random index
        beta <- list()
        for (i in all_indices) {
          beta[[paste0("beta", i)]] <- rnorm(n_features_one, 0, 0.05)
        }

        # Assign corresponding values to betas variables based on assigned_indices_features
        for (i in all_indices) {
          indices <- assigned_indices_features[[i]]
          beta[[i]][indices] <- rnorm(length(indices), 0, 0.05)  # Adjust values dynamically
        }

        pattern_alpha <-"^alpha\\d+"
        matches_alpha <- grepl(pattern_alpha, names(alpha), perl = TRUE)
        n_alpha = length(matches_alpha[matches_alpha == TRUE])

        pattern_beta <- "^beta\\d+"
        matches_beta <- grepl(pattern_beta, names(beta), perl = TRUE)
        n_beta = length(matches_beta[matches_beta == TRUE])

        # Initialize the data list to store the results
        data_list_i <- list()

        # Loop through each alpha and beta combination
        for (i in 1:min(length(alpha), length(beta))) {
          data_i <- alpha[[i]] %*% t(beta[[i]][1:n_features_one])
          data_list_i[[paste0("data.", i)]] <- data_i
        }

        # Combine the results into a single data variable
        data.1 <- Reduce(`+`, data_list_i)

        omic.one <- list()

        # Noise in the first detasets
        sigma <- sigmas_vector[1]

        eps1 <- rnorm(n_samples*n_features_one, 0, sigma) # noise
        omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
        colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
        rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

        gamma <- list()
        gamma <- alpha #list(), assign the values of alpha for the second data

        assigned_indices_features_omic.two <- divide_features_two(n_features_two, num.factor = 1)

        all_indices <- seq_along(assigned_indices_features_omic.two)

        # Create only one vector in list_betas using the random index
        delta <- list()
        for (i in all_indices) {
          delta[[paste0("delta", i)]] <- rnorm(n_features_two, 0, 0.05)
        }

        # Assign corresponding values to betas variables based on assigned_indices_features
        for (i in all_indices) {
          indices <- assigned_indices_features_omic.two[[i]]
          delta[[i]][indices] <- rnorm(length(indices), 10 + 0.05*i, 0.05)  # Adjust values dynamically
        }

        pattern_gamma <-"^alpha\\d+"
        matches_gamma <- grepl(pattern_gamma, names(gamma), perl = TRUE)
        n_gamma = length(matches_gamma[matches_gamma == TRUE])

        pattern_delta <- "^delta\\d+"
        matches_delta <- grepl(pattern_delta, names(delta), perl = TRUE)
        n_delta = length(matches_delta[matches_delta == TRUE])

        # Initialize the data list to store the results
        data_list_j <- list()

        # Loop through each alpha and beta combination
        for (j in 1:min(length(gamma), length(delta))) {
          data_j <- gamma[[j]] %*% t(delta[[j]][1:n_features_two])
          data_list_j[[paste0("data.", j)]] <- data_j
        }

        # Combine the results into a single data variable
        data.2 <- Reduce(`+`, data_list_j)

        omic.two <- list()

        # Noise in the second detasets
        sigma <- sigmas_vector[2]

        eps2 <- rnorm(n_samples*n_features_two, 0, sigma) # noise
        omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
        colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
        rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

        simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

        concatenated_datasets <- list()
        for (i in 1:length(simulated_datasets$object.two)) {
          concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
          concatenated_datasets[[i]] <- concatenated_data
        }
        sim_output <- return(list(concatenated_datasets = concatenated_datasets, divide_samples = divide_samples, assigned_indices_features=assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, list_alphas=alpha, list_gammas=gamma, list_betas=beta, list_deltas=delta))
        return(sim_output)
      }

    }
    else if (type.factor == 'shared'){
      alpha <- list()
      for (i in seq_along(num)) {
        alpha[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
      }

      # Assigning the signal to the selected indices
      for (i in seq_along(alpha)) {
        indices <- assigned_indices_samples
        alpha[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
      }

      # First OMIC data
      select_consecutive_elements_features_from_start <- function(n_features_one) {
        # Ensure the vector has at least 10% elements
        vector = 1:n_features_one

        if (length(vector) < 0.1*n_features_one) {
          stop("The vector must have at least 10% elements.")
        }

        # Determine the number of elements to select (between 10 and the smaller of 90 or the length of the vector)
        num_elements <- sample(0.1*n_features_one:min(0.5*n_features_one, length(vector)), 1)

        # Select the consecutive elements starting from the first element
        selected_elements <- vector[1:num_elements]
        return(selected_elements)
      }

      # FeaturesD signal start_and_end
      assigned_indices_features = select_consecutive_elements_features_from_start(n_features_one=n_features_one)#divide_features_one(n_features_one, num.factor = 1)    # feature indices

      #all_indices <- seq_along(num)

      # Create only one vector in list_betas using the random index
      beta <- list()
      for (i in seq_along(num)) {
        beta[[paste0("beta", i)]] <- rnorm(n_features_one, 0, 0.05)
      }

      # Assign corresponding values to betas variables based on assigned_indices_features
      #all_indices <- seq_along(assigned_indices_features)
      for (i in seq_along(num)) {
        indices <- assigned_indices_features#[[i]]
        beta[[i]][indices] <- rnorm(length(indices), (4.0 + 0.5 * i), 0.05)  # Adjust values dynamically
      }

      pattern_alpha <-"^alpha\\d+"
      matches_alpha <- grepl(pattern_alpha, names(alpha), perl = TRUE)
      n_alpha = length(matches_alpha[matches_alpha == TRUE])

      pattern_beta <- "^beta\\d+"
      matches_beta <- grepl(pattern_beta, names(beta), perl = TRUE)
      n_beta = length(matches_beta[matches_beta == TRUE])

      # Initialize the data list to store the results
      data_list_i <- list()

      # Loop through each alpha and beta combination
      for (i in 1:min(length(alpha), length(beta))) {
        data_i <- alpha[[i]] %*% t(beta[[i]][1:n_features_one])
        data_list_i[[paste0("data.", i)]] <- data_i
      }

      # Combine the results into a single data variable
      data.1 <- Reduce(`+`, data_list_i)

      omic.one <- list()

      # Noise in the first detasets
      sigma <- sigmas_vector[1]

      eps1 <- rnorm(n_samples*n_features_one, 0, sigma) # noise
      omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
      colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
      rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

      # Second OMIC
      # Generate random gamma values based on the max_factors
      gamma <- list()
      gamma <- alpha #list(), assign the values of alpha for the second data

      select_consecutive_elements_features_ending_with_last <- function(n_features_two) {
        # Ensure the vector has at least 10% elements
        vector = 1:n_features_two
        if (length(vector) < 0.1*n_features_two) {
          stop("The vector must have at least 10% elements.")
        }

        # Determine the number of elements to select (between 10 and the smaller of 90 or the length of the vector)
        num_elements <- sample(0.1*n_features_two:min(0.6*n_features_two, length(vector)), 1)

        # Select the consecutive elements ending at the last element
        selected_elements <- vector[(length(vector) - num_elements + 1):(length(vector)+1)]

        return(selected_elements)
      }

      assigned_indices_features_omic.two <- select_consecutive_elements_features_ending_with_last(n_features_two=n_features_two)#divide_features_two(n_features_two, num.factor = 1)

      # all_indices <- seq_along(assigned_indices_features_omic.two)

      # Create only one vector in list_betas using the random index
      delta <- list()
      for (i in seq_along(num)) {
        delta[[paste0("delta", i)]] <- rnorm(n_features_two, 0, 0.05)
      }

      # Assign corresponding values to betas variables based on assigned_indices_features
      for (i in seq_along(num)) {
        indices <- assigned_indices_features_omic.two#[[i]]
        delta[[i]][indices] <- rnorm(length(indices), 6 + 0.05*i, 0.05)  # Adjust values dynamically
      }

      pattern_gamma <-"^alpha\\d+"
      matches_gamma <- grepl(pattern_gamma, names(gamma), perl = TRUE)
      n_gamma = length(matches_gamma[matches_gamma == TRUE])

      pattern_delta <- "^delta\\d+"
      matches_delta <- grepl(pattern_delta, names(delta), perl = TRUE)
      n_delta = length(matches_delta[matches_delta == TRUE])

      # Initialize the data list to store the results
      data_list_j <- list()

      # Loop through each alpha and beta combination
      for (j in 1:min(length(gamma), length(delta))) {
        data_j <- gamma[[j]] %*% t(delta[[j]][1:n_features_two])
        data_list_j[[paste0("data.", j)]] <- data_j
      }

      # Combine the results into a single data variable
      data.2 <- Reduce(`+`, data_list_j)

      omic.two <- list()

      # Noise in the second detasets
      sigma <- sigmas_vector[2]

      eps2 <- rnorm(n_samples*n_features_two, 0, sigma) # noise
      omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
      colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
      rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

      simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

      concatenated_datasets <- list()
      for (i in 1:length(simulated_datasets$object.two)) {
        concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
        concatenated_datasets[[i]] <- concatenated_data
      }
      sim_output <- return(list(concatenated_datasets = concatenated_datasets, divide_samples = divide_samples, assigned_indices_features=assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, list_alphas=alpha, list_gammas=gamma, list_betas=beta, list_deltas=delta))
      return(sim_output)

    }else if (type.factor==""){
      selected_factor <- c('unique', 'shared')
      selected.factor <- sample (selected_factor, 1)

      if (selected.factor == 'unique'){
        omic <- c('omic.one', 'omic.two')
        selected.omic <- sample(omic, 1)

        if (selected.omic == 'omic.one'){
          alpha <- list()
          for (i in seq_along(assigned_indices_samples)) {
            alpha[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
          }

          # Assigning the signal to the selected indices
          for (i in seq_along(alpha)) {
            indices <- assigned_indices_samples#[[i]]
            alpha[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
          }

          # # First OMIC data
          # FeaturesD signal start_and_end
          assigned_indices_features = divide_features_one(n_features_one, num.factor = 1)    # feature indices

          all_indices <- seq_along(assigned_indices_features)

          # Create only one vector in list_betas using the random index
          beta <- list()
          for (i in all_indices) {
            beta[[paste0("beta", i)]] <- rnorm(n_features_one, 0, 0.05)
          }

          # Assign corresponding values to betas variables based on assigned_indices_features
          for (i in all_indices) {
            indices <- assigned_indices_features[[i]]
            beta[[i]][indices] <- rnorm(length(indices), (4.0 + 0.5 * i), 0.05)  # Adjust values dynamically
          }

          pattern_alpha <-"^alpha\\d+"
          matches_alpha <- grepl(pattern_alpha, names(alpha), perl = TRUE)
          n_alpha = length(matches_alpha[matches_alpha == TRUE])

          pattern_beta <- "^beta\\d+"
          matches_beta <- grepl(pattern_beta, names(beta), perl = TRUE)
          n_beta = length(matches_beta[matches_beta == TRUE])

          # Initialize the data list to store the results
          data_list_i <- list()

          # Loop through each alpha and beta combination
          for (i in 1:min(length(alpha), length(beta))) {
            data_i <- alpha[[i]] %*% t(beta[[i]][1:n_features_one])
            data_list_i[[paste0("data.", i)]] <- data_i
          }

          # Combine the results into a single data variable
          data.1 <- Reduce(`+`, data_list_i)

          omic.one <- list()

          # Noise in the first detasets
          sigma <- sigmas_vector[1]

          eps1 <- rnorm(n_samples*n_features_one, 0, sigma) # noise
          omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
          colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
          rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

          # Second OMIC
          # Generate random gamma values based on the max_factors
          gamma <- list()
          gamma <- alpha #list(), assign the values of alpha for the second data

          assigned_indices_features_omic.two <- divide_features_two(n_features_two, num.factor = 1)

          all_indices <- seq_along(assigned_indices_features_omic.two)

          # Create only one vector in list_betas using the random index
          delta <- list()
          for (i in all_indices) {
            delta[[paste0("delta", i)]] <- rnorm(n_features_two, 0, 0.05)
          }

          # Assign corresponding values to betas variables based on assigned_indices_features
          for (i in all_indices) {
            indices <- assigned_indices_features_omic.two[[i]]
            delta[[i]][indices] <- rnorm(length(indices), 0, 0.05)  # Adjust values dynamically
          }

          pattern_gamma <-"^alpha\\d+"
          matches_gamma <- grepl(pattern_gamma, names(gamma), perl = TRUE)
          n_gamma = length(matches_gamma[matches_gamma == TRUE])

          pattern_delta <- "^delta\\d+"
          matches_delta <- grepl(pattern_delta, names(delta), perl = TRUE)
          n_delta = length(matches_delta[matches_delta == TRUE])

          # Initialize the data list to store the results
          data_list_j <- list()

          # Loop through each alpha and beta combination
          for (j in 1:min(length(gamma), length(delta))) {
            data_j <- gamma[[j]] %*% t(delta[[j]][1:n_features_two])
            data_list_j[[paste0("data.", j)]] <- data_j
          }

          # Combine the results into a single data variable
          data.2 <- Reduce(`+`, data_list_j)

          omic.two <- list()

          # Noise in the second detasets
          sigma <- sigmas_vector[2]

          eps2 <- rnorm(n_samples*n_features_two, 0, sigma) # noise
          omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
          colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
          rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

          phi <- ceiling(n_features_one)

          simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

          concatenated_datasets <- list()
          for (i in 1:length(simulated_datasets$object.two)) {
            concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
            concatenated_datasets[[i]] <- concatenated_data
          }
          sim_output <- return(list(concatenated_datasets = concatenated_datasets, divide_samples = divide_samples, assigned_indices_features=assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, list_alphas=alpha, list_gammas=gamma, list_betas=beta, list_deltas=delta))
          return(sim_output)

        }else if (selected.omic == 'omic.two') {
          alpha <- list()
          for (i in seq_along(assigned_indices_samples)) {
            alpha[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
          }

          # Assigning the signal to the selected indices
          for (i in seq_along(alpha)) {
            indices <- assigned_indices_samples#[[i]]
            alpha[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
          }


          assigned_indices_features = divide_features_one(n_features_one, num.factor = 1)

          all_indices <- seq_along(assigned_indices_features)

          # Create only one vector in list_betas using the random index
          beta <- list()
          for (i in all_indices) {
            beta[[paste0("beta", i)]] <- rnorm(n_features_one, 0, 0.05)
          }

          # Assign corresponding values to betas variables based on assigned_indices_features
          for (i in all_indices) {
            indices <- assigned_indices_features[[i]]
            beta[[i]][indices] <- rnorm(length(indices), 0, 0.05)  # Adjust values dynamically
          }

          pattern_alpha <-"^alpha\\d+"
          matches_alpha <- grepl(pattern_alpha, names(alpha), perl = TRUE)
          n_alpha = length(matches_alpha[matches_alpha == TRUE])

          pattern_beta <- "^beta\\d+"
          matches_beta <- grepl(pattern_beta, names(beta), perl = TRUE)
          n_beta = length(matches_beta[matches_beta == TRUE])

          # Initialize the data list to store the results
          data_list_i <- list()

          # Loop through each alpha and beta combination
          for (i in 1:min(length(alpha), length(beta))) {
            data_i <- alpha[[i]] %*% t(beta[[i]][1:n_features_one])
            data_list_i[[paste0("data.", i)]] <- data_i
          }

          # Combine the results into a single data variable
          data.1 <- Reduce(`+`, data_list_i)

          omic.one <- list()

          # Noise in the first detasets
          sigma <- sigmas_vector[1]

          eps1 <- rnorm(n_samples*n_features_one, 0, sigma) # noise
          omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
          colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
          rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

          gamma <- list()
          gamma <- alpha #list(), assign the values of alpha for the second data

          assigned_indices_features_omic.two <- divide_features_two(n_features_two, num.factor = 1)

          #all_indices <- seq_along(assigned_indices_features_omic.two)

          # Create only one vector in list_betas using the random index
          delta <- list()
          for (i in seq_along(num)) {
            delta[[paste0("delta", i)]] <- rnorm(n_features_two, 0, 0.05)
          }

          # Assign corresponding values to betas variables based on assigned_indices_features
          #for (i in seq_along(num)) {
          #  indices <- assigned_indices_features_omic.two#[[i]]
          #  delta[[i]][indices] <- rnorm(length(indices), 7 + 0.05*i, 0.05)  # Adjust values dynamically
          #}
          # Assign corresponding values to betas variables based on assigned_indices_features
          for (i in seq_along(num)) {
            # Extract indices correctly; if indices is a list, unlist it to make it a vector
            indices <- unlist(assigned_indices_features_omic.two)  # Use [[i]] if you need a specific subset

            # Ensure indices are not NA and are numeric/integer
            indices <- indices[!is.na(indices)]

            # Check if indices are numeric and not empty
            if (length(indices) > 0 && is.numeric(indices)) {
              # Assign random values to delta[[i]] at the specified indices
              delta[[i]][indices] <- rnorm(length(indices), 7 + 0.05 * i, 0.05)  # Adjust values dynamically
            } else {
              warning(paste("Invalid or empty indices for iteration:", i))
            }
          }

          pattern_gamma <-"^alpha\\d+"
          matches_gamma <- grepl(pattern_gamma, names(gamma), perl = TRUE)
          n_gamma = length(matches_gamma[matches_gamma == TRUE])

          pattern_delta <- "^delta\\d+"
          matches_delta <- grepl(pattern_delta, names(delta), perl = TRUE)
          n_delta = length(matches_delta[matches_delta == TRUE])

          # Initialize the data list to store the results
          data_list_j <- list()

          # Loop through each alpha and beta combination
          for (j in 1:min(length(gamma), length(delta))) {
            data_j <- gamma[[j]] %*% t(delta[[j]][1:n_features_two])
            data_list_j[[paste0("data.", j)]] <- data_j
          }

          # Combine the results into a single data variable
          data.2 <- Reduce(`+`, data_list_j)

          omic.two <- list()

          # Noise in the second detasets
          sigma <- sigmas_vector[2]

          eps2 <- rnorm(n_samples*n_features_two, 0, sigma) # noise
          omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
          colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
          rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

          simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

          concatenated_datasets <- list()
          for (i in 1:length(simulated_datasets$object.two)) {
            concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
            concatenated_datasets[[i]] <- concatenated_data
          }
          sim_output <- return(list(concatenated_datasets = concatenated_datasets, divide_samples = divide_samples, assigned_indices_features=assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, list_alphas=alpha, list_gammas=gamma, list_betas=beta, list_deltas=delta))
          return(sim_output)
        }
      }else if (selected.factor == 'shared'){
        alpha <- list()
        for (i in seq_along(num)) {
          alpha[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
        }

        # Assigning the signal to the selected indices
        for (i in seq_along(alpha)) {
          indices <- assigned_indices_samples
          alpha[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
        }

        # First OMIC data
        select_consecutive_elements_features_from_start <- function(n_features_one) {
          # Ensure the vector has at least 10% elements
          vector = 1:n_features_one

          if (length(vector) < 0.1*n_features_one) {
            stop("The vector must have at least 10% elements.")
          }

          # Determine the number of elements to select (between 10 and the smaller of 90 or the length of the vector)
          num_elements <- sample(0.1*n_features_one:min(0.5*n_features_one, length(vector)), 1)

          # Select the consecutive elements starting from the first element
          selected_elements <- vector[1:num_elements]
          return(selected_elements)
        }

        # FeaturesD signal start_and_end
        assigned_indices_features = select_consecutive_elements_features_from_start(n_features_one=n_features_one)#divide_features_one(n_features_one, num.factor = 1)    # feature indices

        #all_indices <- seq_along(num)

        # Create only one vector in list_betas using the random index
        beta <- list()
        for (i in seq_along(num)) {
          beta[[paste0("beta", i)]] <- rnorm(n_features_one, 0, 0.05)
        }

        # Assign corresponding values to betas variables based on assigned_indices_features
        #all_indices <- seq_along(assigned_indices_features)
        for (i in seq_along(num)) {
          indices <- assigned_indices_features#[[i]]
          beta[[i]][indices] <- rnorm(length(indices), (4.0 + 0.5 * i), 0.05)  # Adjust values dynamically
        }

        pattern_alpha <-"^alpha\\d+"
        matches_alpha <- grepl(pattern_alpha, names(alpha), perl = TRUE)
        n_alpha = length(matches_alpha[matches_alpha == TRUE])

        pattern_beta <- "^beta\\d+"
        matches_beta <- grepl(pattern_beta, names(beta), perl = TRUE)
        n_beta = length(matches_beta[matches_beta == TRUE])

        # Initialize the data list to store the results
        data_list_i <- list()

        # Loop through each alpha and beta combination
        for (i in 1:min(length(alpha), length(beta))) {
          data_i <- alpha[[i]] %*% t(beta[[i]][1:n_features_one])
          data_list_i[[paste0("data.", i)]] <- data_i
        }

        # Combine the results into a single data variable
        data.1 <- Reduce(`+`, data_list_i)

        omic.one <- list()

        # Noise in the first detasets
        sigma <- sigmas_vector[1]

        eps1 <- rnorm(n_samples*n_features_one, 0, sigma) # noise
        omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
        colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
        rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

        # Second OMIC
        # Generate random gamma values based on the max_factors
        gamma <- list()
        gamma <- alpha #list(), assign the values of alpha for the second data

        select_consecutive_elements_features_ending_with_last <- function(n_features_two) {
          # Ensure the vector has at least 10% elements
          vector = 1:n_features_two
          if (length(vector) < 0.1*n_features_two) {
            stop("The vector must have at least 10% elements.")
          }

          # Determine the number of elements to select (between 10 and the smaller of 90 or the length of the vector)
          num_elements <- sample(0.1*n_features_two:min(0.6*n_features_two, length(vector)), 1)

          # Select the consecutive elements ending at the last element
          selected_elements <- vector[(length(vector) - num_elements + 1):(length(vector)+1)]

          return(selected_elements)
        }

        assigned_indices_features_omic.two <- select_consecutive_elements_features_ending_with_last(n_features_two=n_features_two)#divide_features_two(n_features_two, num.factor = 1)

        all_indices <- seq_along(assigned_indices_features_omic.two)

        # Create only one vector in list_betas using the random index
        delta <- list()
        for (i in seq_along(num)) {
          delta[[paste0("delta", i)]] <- rnorm(n_features_two, 0, 0.05)
        }

        # Assign corresponding values to betas variables based on assigned_indices_features
        # Loop through each element in num
        for (i in seq_along(num)) {
          # Extract indices, ensuring no NA values are present
          indices <- assigned_indices_features_omic.two #[[i]] # If you meant to subset using [[i]], ensure it's correct

          # Remove NAs from indices to avoid assignment errors
          indices <- indices[!is.na(indices)]

          # Check if indices is not empty after removing NAs
          if (length(indices) > 0) {
            # Assign values to delta[[i]] at the specified indices
            delta[[i]][indices] <- rnorm(length(indices), 8 + 0.05 * i, 0.05)  # Adjust values dynamically
          } else {
            warning(paste("No valid indices for iteration:", i)) # Optional warning if indices are empty
          }
        }

        pattern_gamma <-"^alpha\\d+"
        matches_gamma <- grepl(pattern_gamma, names(gamma), perl = TRUE)
        n_gamma = length(matches_gamma[matches_gamma == TRUE])

        pattern_delta <- "^delta\\d+"
        matches_delta <- grepl(pattern_delta, names(delta), perl = TRUE)
        n_delta = length(matches_delta[matches_delta == TRUE])

        # Initialize the data list to store the results
        data_list_j <- list()

        # Loop through each alpha and beta combination
        for (j in 1:min(length(gamma), length(delta))) {
          data_j <- gamma[[j]] %*% t(delta[[j]][1:n_features_two])
          data_list_j[[paste0("data.", j)]] <- data_j
        }

        # Combine the results into a single data variable
        data.2 <- Reduce(`+`, data_list_j)

        omic.two <- list()

        # Noise in the second detasets
        sigma <- sigmas_vector[2]

        eps2 <- rnorm(n_samples*n_features_two, 0, sigma) # noise
        omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
        colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
        rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

        simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

        concatenated_datasets <- list()
        for (i in 1:length(simulated_datasets$object.two)) {
          concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
          concatenated_datasets[[i]] <- concatenated_data
        }
        sim_output <- return(list(concatenated_datasets = concatenated_datasets, divide_samples = divide_samples, assigned_indices_features=assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, list_alphas=alpha, list_gammas=gamma, list_betas=beta, list_deltas=delta))
        return(sim_output)
      }
    }
  }else if(num.factor == 'multiple'){
    if(is.null(advanced_dist) || advanced_dist == '' || advanced_dist == 'mixed'){
      # Code to execute if advanced_dist is either an empty string NULL, "" or "mixed"

      # features
      n_features_one <- vector_features[1]
      n_features_two <- vector_features[2]

      # Noise in the first detasets
      sigmas <- sigmas_vector[1]

      #num_factor = num_factors
      omic.one <- list()
      omic.two <- list()

      #assigned_indices_samples <- divide_samples(n_samples = n_samples, num=n_factors, min_size=5)
      if (n_samples <= 20) {
        min_size <- 2
        num = n_factors
        check <- num * min_size
        # Ensure that the total minimum size does not fall below 75% of n_samples
        if (check > n_samples * 0.75) {
          stop("The minimum segment size is too small given the constraints.")  # Replace with the appropriate error message
        } else {
          assigned_indices_samples <- divide_samples(n_samples = n_samples, num = num, min_size = min_size)  # min_size is set within divide_samples
        }
      } else if (n_samples > 20) {
        min_size <- 5
        num = n_factors
        check <- num * min_size
        # Ensure that the total minimum size does not fall below 75% of n_samples
        if (check > n_samples * 0.75) {
          stop("The minimum segment size is too small given the constraints.")  # Replace with the appropriate error message
        } else {
          assigned_indices_samples <- divide_samples(n_samples = n_samples, num = num, min_size = min_size)  # min_size is set within divide_samples
        }
      }

      #print(assigned_indices_samples)

      # sample indices
      # Check the max number of factors in assigned_indices_samples
      max_factors <- length(assigned_indices_samples)# Same as n_factors

      # Generate random alpha values based on the max_factors
      list_alphas <- list()
      for (i in 1:max_factors) {
        list_alphas[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
      }

      # Assign corresponding values to alpha variables based on assigned_indices_samples
      for (i in seq_along(assigned_indices_samples)) {
        indices <- assigned_indices_samples[[i]]
        list_alphas[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
      }
      list_gammas <- list_alphas
      # Add the separated list of the indices selected for omics of multiple features
      # Initial vector
      vector <- c(1:max_factors)

      # Specify the number of elements to sample for the shared vector
      num_shared <- ceiling(runif(1, min = 1, max = length(vector) - 1))

      # Select the 'shared' elements
      if (num_shared > 0) {
        shared <- sample(vector, num_shared)
        # Remove the selected shared elements from the vector
        remain_vector <- vector[!vector %in% shared]
      }

      # Shuffle the remaining vector
      if(length(remain_vector == 1)){
        shuffled_remain = remain_vector
      } else{
        shuffled_remain <- sample(remain_vector)
      }
      # Split the shuffled elements equally or nearly equally between omic_one_unique and omic_two_unique
      num_elements <- length(shuffled_remain)
      split_point <- sample(0:num_elements, 1) # Randomly select a split point

      omic_one_unique <- shuffled_remain[0:split_point] # First part goes to omic_one_unique
      omic_two_unique <- shuffled_remain[!shuffled_remain %in% omic_one_unique] # Second part goes to omic_two_unique

      # Print the results
      list(
        shared = shared,
        omic_one_unique = omic_one_unique,
        omic_two_unique = omic_two_unique
      )

      # End of indices assignment

      # Assignment of factors
      list_omic_one_factors = c(shared, omic_one_unique)
      list_omic_two_factors = c(shared, omic_two_unique)
      #all_indices = c(shared, omic_one_unique, omic_two_unique)

      factor_xtics <- list(
        shared = shared,
        omic_one_unique = omic_one_unique,
        omic_two_unique = omic_two_unique
      )

      ### * * * * * * * * * * * * * * * * * * * * ** * * * * First OMIC data * * * * * * * * * * * * * * * * * * * * * * * * * *#

      # FeaturesD signal start_and_end
      assigned_indices_features = feature_selection_one(n_features_one = n_features_one, num.factor = "multiple", no_factor = length(list_omic_one_factors))

      # Generate random beta values based on the max_factors
      list_betas <- list()

      # Shuffle the list
      (shuffled_assigned_indices_features <- assigned_indices_features)# sample(assigned_indices_features)) # Do not reshuffle first

      # Create only one vector in list_betas using the random index
      for (i in seq_along(assigned_indices_features)) {
        list_betas[[i]] <- rnorm(n_features_one, 0, 0.05)
      }

      # Loop through the ordered indices in list_omic_one_factors
      for (i in seq_along(assigned_indices_features)) {
        indices_ns <- assigned_indices_features[[i]]  # Get the corresponding indices from assigned_indices_features
        if (length(indices_ns) > 0) {
          # Initialize list_betas[[i]] if it's not already
          if (is.null(list_betas[[i]])) {
            list_betas[[i]] <- numeric(max(indices_ns))
          }
          # Assign values to the appropriate indices in list_betas
          list_betas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (5.0 + 0.4 * i), sd = 0.05)
        }
      }

      # Renaming list_betas to desired factor number
      for (i in seq_along(list_omic_one_factors)) {
        names(list_betas)[i] <- paste0("beta", list_omic_one_factors[[i]])#list_omic_one_factors[[i]]
      }

      pattern_alpha <-"^alpha\\d+"
      matches_alpha <- grepl(pattern_alpha, names(list_alphas), perl = TRUE)
      n_alpha = length(matches_alpha[matches_alpha == TRUE])

      pattern_beta <- "^beta\\d+"
      matches_beta <- grepl(pattern_beta, names(list_betas), perl = TRUE)
      n_beta = length(matches_beta[matches_beta == TRUE])

      # Initialize the data list to store the results
      data_list_i <- list()

      # Sort the list by names in ascending order
      list_betas <- list_betas[order(names(list_betas))]
      list_alphas <- list_alphas[order(names(list_alphas))]

      # Extract the last numeric values from the names of the vectors
      alphas_names <- as.numeric(gsub("[^0-9]", "", names(list_alphas)))
      betas_names <- as.numeric(gsub("[^0-9]", "", names(list_betas)))

      # Find common numbers between the last values of the vector names
      common_names <- intersect(alphas_names, betas_names)

      # Filter the vectors in each list based on the common names
      list_alphas <- list_alphas[paste0("alpha", common_names)]
      list_betas <- list_betas[paste0("beta", common_names)]

      # Loop through each alpha and beta combination
      for (i in 1:min(length(list_alphas), length(list_betas))) {
        data_i <- list_alphas[[i]] %*% t(list_betas[[i]][1:n_features_one])
        data_list_i[[paste0("data.", i)]] <- data_i
      }

      # Combine the results into a single data variable
      data.1 <- Reduce(`+`, data_list_i)

      eps1 <- rnorm(n_samples * n_features_one, 0, sigmas) # noise
      omic.one <- list()
      omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
      colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
      rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

      #dataset <- omic.one[[1]]
      #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")

      ### * * * * * * * * * * * * * * * * * * * * ** * * * * Second OMIC data * * * * * * * * * * * * * * * * * * * * * * * * * *

      # Noise in the first detasets
      sigmas <- sigmas_vector[2]

      # Generate random gamma values based on the max_factors
      # Replace 'alpha' with 'gamma' in the names of list_gammas
      names(list_gammas) <- gsub("alpha", "gamma", names(list_gammas))

      # Extract the numeric parts from the names of list_gammas
      numeric_part <- as.numeric(gsub("\\D", "", names(list_gammas)))

      # Retain elements where the numeric part of the name matches values in list_omic_two_factors
      list_gammas <- list_gammas[numeric_part %in% list_omic_two_factors]

      # Features
      assigned_indices_features_omic.two = feature_selection_two(n_features_two = n_features_two, num.factor = "multiple", no_factor = length(list_omic_two_factors))

      # Empty list
      list_deltas <- list()

      # Create only one vector in list_deltas using the random index
      for (i in seq_along(assigned_indices_features_omic.two)) {
        list_deltas[[i]] <- rnorm(n_features_two, 0, 0.05)
      }

      # # Loop through the ordered indices in list_omic_one_factors
      for (i in seq_along(assigned_indices_features_omic.two)) {
        indices_ns <- assigned_indices_features_omic.two[[i]]  # Get the corresponding indices from assigned_indices_features_omic.two
        # Remove NA values from indices_ns
        indices_ns <- indices_ns[!is.na(indices_ns)]
        if (length(indices_ns) > 0) {
          # Initialize list_deltas[[i]] if it's not already
          if (is.null(list_deltas[[i]])) {
            list_deltas[[i]] <- numeric(max(indices_ns, na.rm = TRUE))  # Use na.rm = TRUE to handle potential NA
          }
          # Ensure indices are within the range of list_deltas[[i]]
          if (max(indices_ns) > length(list_deltas[[i]])) {
            length(list_deltas[[i]]) <- max(indices_ns)
          }
          # Assign values to the appropriate indices in list_deltas
          list_deltas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (5.0 + 0.4 * i), sd = 0.05)
        }
      }

      # Renaming list_deltas to desired factor number
      for (i in seq_along(list_omic_two_factors)) {
        names(list_deltas)[i] <- paste0("delta", list_omic_two_factors[[i]])#list_omic_one_factors[[i]]
      }

      pattern_gamma <-"^gamma\\d+"
      matches_gamma <- grepl(pattern_gamma, names(list_gammas), perl = TRUE)
      n_gamma = length(matches_alpha[matches_gamma == TRUE])

      pattern_delta <- "^delta\\d+"
      matches_delta <- grepl(pattern_delta, names(list_deltas), perl = TRUE)
      n_delta = length(matches_delta[matches_delta == TRUE])

      # Initialize the data list to store the results
      data_list_j <- list()

      # Sort the list by names in ascending order
      list_gammas <- list_gammas[order(names(list_gammas))]
      list_deltas <- list_deltas[order(names(list_deltas))]

      # Extract the last numeric values from the names of the vectors
      gammas_names <- as.numeric(gsub("[^0-9]", "", names(list_gammas)))
      deltas_names <- as.numeric(gsub("[^0-9]", "", names(list_deltas)))

      # Find common numbers between the last values of the vector names
      common_names <- intersect(gammas_names, deltas_names)

      # Filter the vectors in each list based on the common names
      list_gammas <- list_gammas[paste0("gamma", common_names)]
      list_deltas <- list_deltas[paste0("delta", common_names)]

      # Loop through each alpha and beta combination
      for (j in 1:min(length(list_gammas), length(list_deltas))) {
        data_j <- list_gammas[[j]] %*% t(list_deltas[[j]][1:n_features_two])
        data_list_j[[paste0("data.", j)]] <- data_j
      }

      # Combine the results into a single data variable
      data.2 <- Reduce(`+`, data_list_j)

      eps2 <- rnorm(n_samples * n_features_two, 0,sigmas)
      omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
      colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
      rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

      #dataset <- data.2#dataset <- omic.two[[1]]
      #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")

      simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

      concatenated_datasets <- list()
      for (i in 1:length(simulated_datasets$object.two)) {
        concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
        concatenated_datasets[[i]] <- concatenated_data
      }

      # sim_output <- list(concatenated_datasets = concatenated_datasets, list_alphas = list_alphas, list_betas = list_betas, list_gammas = list_gammas, list_deltas = list_deltas)#simulated_datasets=simulated_datasets)
      sim_output <- list(concatenated_datasets=concatenated_datasets, omic.one=omic.one, omic.two=omic.two, assigned_indices_samples = assigned_indices_samples, assigned_indices_features = assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, factor_xtics=factor_xtics, list_alphas=list_alphas, list_gammas=list_gammas, list_betas=list_betas, list_deltas=list_deltas)#simulated_datasets=simulated_datasets)

      return(sim_output)

    }else if(advanced_dist == 'omic.one'){
      # Code to execute if advanced_dist is "omic.one"

      # features
      n_features_one <- vector_features[1]
      n_features_two <- vector_features[2]

      # Noise in the first detasets
      sigmas <- sigmas_vector[1]

      #num_factor = num_factors
      omic.one <- list()
      omic.two <- list()

      assigned_indices_samples <- divide_samples(n_samples = n_samples, num=n_factors, min_size=10)
      #print(assigned_indices_samples)

      # sample indices
      # Check the max number of factors in assigned_indices_samples
      max_factors <- length(assigned_indices_samples)# Same as n_factors

      # Generate random alpha values based on the max_factors
      list_alphas <- list()
      for (i in 1:max_factors) {
        list_alphas[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
      }

      # Assign corresponding values to alpha variables based on assigned_indices_samples
      for (i in seq_along(assigned_indices_samples)) {
        indices <- assigned_indices_samples[[i]]
        list_alphas[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
      }
      list_gammas <- list_alphas
      # Add the separated list of the indices selected for omics of multiple features
      # Initial vector
      vector <- c(1:max_factors)

      # Specify the number of elements to sample for the shared vector
      num_shared <- 0#ceiling(runif(1, min = 1, max = length(vector) - 1))

      # Select the 'shared' elements
      if (num_shared > 0) {
        shared <- sample(vector, num_shared)
        # Remove the selected shared elements from the vector
        remain_vector <- vector[!vector %in% shared]
      }else{
        shared <- 0
        remain_vector <- vector
      }

      # Shuffle the remaining vector
      if(length(remain_vector == 1)){
        shuffled_remain = remain_vector
      } else{
        shuffled_remain <- sample(remain_vector)
      }
      # Split the shuffled elements equally or nearly equally between omic_one_unique and omic_two_unique
      num_elements <- length(shuffled_remain)
      split_point <- sample(0:num_elements, 1) # Randomly select a split point

      omic_one_unique <- remain_vector #shuffled_remain[0:split_point] # First part goes to omic_one_unique
      omic_two_unique <- 0#shuffled_remain[!shuffled_remain %in% omic_one_unique] # Second part goes to omic_two_unique

      # Print the results
      list(
        shared = shared,
        omic_one_unique = omic_one_unique,
        omic_two_unique = omic_two_unique
      )

      # End of indices assignment

      # Assignment of factors
      list_omic_one_factors = c(shared, omic_one_unique)
      list_omic_two_factors = c(shared, omic_two_unique)
      #all_indices = c(shared, omic_one_unique, omic_two_unique)

      factor_xtics <- list(
        shared = shared,
        omic_one_unique = omic_one_unique,
        omic_two_unique = omic_two_unique
      )

      ### * * * * * * * * * * * * * * * * * * * * ** * * * * First OMIC data * * * * * * * * * * * * * * * * * * * * * * * * * *#

      # FeaturesD signal start_and_end
      list_omic_factors <- seq_along(1:n_factors)
      assigned_indices_features = feature_selection_one(n_features_one = n_features_one, num.factor = "multiple", no_factor = length(list_omic_factors))

      # Generate random beta values based on the max_factors
      list_betas <- list()

      # Shuffle the list
      (shuffled_assigned_indices_features <- sample(assigned_indices_features)) # assigned_indices_features)# Do not reshuffle first

      # Create only one vector in list_betas using the random index
      for (i in seq_along(shuffled_assigned_indices_features)) {
        list_betas[[i]] <- rnorm(n_features_one, 0, 0.05)
      }

      # Loop through the ordered indices in list_omic_one_factors
      for (i in seq_along(shuffled_assigned_indices_features)) {
        indices_ns <- shuffled_assigned_indices_features[[i]]  # Get the corresponding indices from assigned_indices_features
        if (length(indices_ns) > 0) {
          # Initialize list_betas[[i]] if it's not already
          if (is.null(list_betas[[i]])) {
            list_betas[[i]] <- numeric(max(indices_ns))
          }
          # Assign values to the appropriate indices in list_betas
          list_betas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (5.0 + 0.4 * i), sd = 0.05)
        }
      }

      # Renaming list_betas to desired factor number
      for (i in seq_along(list_omic_factors)) {
        names(list_betas)[i] <- paste0("beta", list_omic_factors[[i]])#list_omic_one_factors[[i]]
      }

      pattern_alpha <-"^alpha\\d+"
      matches_alpha <- grepl(pattern_alpha, names(list_alphas), perl = TRUE)
      n_alpha = length(matches_alpha[matches_alpha == TRUE])

      pattern_beta <- "^beta\\d+"
      matches_beta <- grepl(pattern_beta, names(list_betas), perl = TRUE)
      n_beta = length(matches_beta[matches_beta == TRUE])

      # Initialize the data list to store the results
      data_list_i <- list()

      # Sort the list by names in ascending order
      list_betas <- list_betas[order(names(list_betas))]
      list_alphas <- list_alphas[order(names(list_alphas))]

      # Extract the last numeric values from the names of the vectors
      alphas_names <- as.numeric(gsub("[^0-9]", "", names(list_alphas)))
      betas_names <- as.numeric(gsub("[^0-9]", "", names(list_betas)))

      # Find common numbers between the last values of the vector names
      common_names <- intersect(alphas_names, betas_names)

      # Filter the vectors in each list based on the common names
      list_alphas <- list_alphas[paste0("alpha", common_names)]
      list_betas <- list_betas[paste0("beta", common_names)]

      # Loop through each alpha and beta combination
      for (i in 1:min(length(list_alphas), length(list_betas))) {
        data_i <- list_alphas[[i]] %*% t(list_betas[[i]][1:n_features_one])
        data_list_i[[paste0("data.", i)]] <- data_i
      }

      # Combine the results into a single data variable
      data.1 <- Reduce(`+`, data_list_i)

      eps1 <- rnorm(n_samples * n_features_one, 0, sigmas) # noise
      omic.one <- list()
      omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
      colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
      rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

      #dataset <- omic.one[[1]]
      #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")

      ### * * * * * * * * * * * * * * * * * * * * ** * * * * Second OMIC data * * * * * * * * * * * * * * * * * * * * * * * * * *

      # Noise in the first detasets
      sigmas <- sigmas_vector[2]

      # Generate random gamma values based on the max_factors
      list_gammas <- list_alphas
      # Replace 'alpha' with 'gamma' in the names of list_gammas
      names(list_gammas) <- gsub("alpha", "gamma", names(list_gammas))

      # Extract the numeric parts from the names of list_gammas
      numeric_part <- as.numeric(gsub("\\D", "", names(list_gammas)))

      # Retain elements where the numeric part of the name matches values in list_omic_two_factors
      list_omic_factors <- seq_along(1:n_factors)
      list_gammas <- list_gammas[numeric_part %in% list_omic_factors]

      # Features
      assigned_indices_features_omic.two = feature_selection_two(n_features_two=n_features_two, num.factor="multiple", no_factor=n_factors)

      # Empty list
      list_deltas <- list()

      # Create only one vector in list_deltas using the random index
      for (i in seq_along(assigned_indices_features_omic.two)) {
        list_deltas[[i]] <- rnorm(n_features_two, 0, 0.05)
      }

      # Loop through the ordered indices in list_omic_one_factors
      for (i in seq_along(assigned_indices_features_omic.two)) {
        indices_ns <- assigned_indices_features_omic.two[[i]]  # Get the corresponding indices from assigned_indices_features
        if (length(indices_ns) > 0) {
          # Initialize list_deltas[[i]] if it's not already
          if (is.null(list_deltas[[i]])) {
            list_deltas[[i]] <- numeric(max(indices_ns))
          }
          # Assign values to the appropriate indices in list_deltas
          list_deltas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = 0, sd = 0.05)
        }
      }

      # Renaming list_deltas to desired factor number
      for (i in seq_along(assigned_indices_features_omic.two)) {
        list_omic_factors = seq_along(assigned_indices_features_omic.two)
        names(list_deltas)[i] <- paste0("delta", list_omic_factors[[i]])#list_omic_one_factors[[i]]
      }

      pattern_gamma <-"^gamma\\d+"
      matches_gamma <- grepl(pattern_gamma, names(list_gammas), perl = TRUE)
      n_gamma = length(matches_alpha[matches_gamma == TRUE])

      pattern_delta <- "^delta\\d+"
      matches_delta <- grepl(pattern_delta, names(list_deltas), perl = TRUE)
      n_delta = length(matches_delta[matches_delta == TRUE])

      # Initialize the data list to store the results
      data_list_j <- list()

      # Sort the list by names in ascending order
      list_gammas <- list_gammas[order(names(list_gammas))]
      list_deltas <- list_deltas[order(names(list_deltas))]

      # Extract the last numeric values from the names of the vectors
      gammas_names <- as.numeric(gsub("[^0-9]", "", names(list_gammas)))
      deltas_names <- as.numeric(gsub("[^0-9]", "", names(list_deltas)))

      # Find common numbers between the last values of the vector names
      common_names <- intersect(gammas_names, deltas_names)

      # Filter the vectors in each list based on the common names
      list_gammas <- list_gammas[paste0("gamma", common_names)]
      list_deltas <- list_deltas[paste0("delta", common_names)]

      # Loop through each alpha and beta combination
      for (j in 1:min(length(list_gammas), length(list_deltas))) {
        data_j <- list_gammas[[j]] %*% t(list_deltas[[j]][1:n_features_two])
        data_list_j[[paste0("data.", j)]] <- data_j
      }

      # Combine the results into a single data variable
      data.2 <- Reduce(`+`, data_list_j)
      omic.two<- list()
      eps2 <- rnorm(n_samples * n_features_two, 0,sigmas)
      omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
      colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
      rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

      #dataset <- data.2#dataset <- omic.two[[1]]
      #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")

      simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

      concatenated_datasets <- list()
      for (i in 1:length(simulated_datasets$object.two)) {
        concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
        concatenated_datasets[[i]] <- concatenated_data
      }

      # sim_output <- list(concatenated_datasets = concatenated_datasets, list_alphas = list_alphas, list_betas = list_betas, list_gammas = list_gammas, list_deltas = list_deltas)#simulated_datasets=simulated_datasets)
      sim_output <- list(concatenated_datasets=concatenated_datasets, omic.one=omic.one, omic.two=omic.two, assigned_indices_samples = assigned_indices_samples, assigned_indices_features = assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, factor_xtics=factor_xtics, list_alphas=list_alphas, list_gammas=list_gammas, list_betas=list_betas, list_deltas=list_deltas)#simulated_datasets=simulated_datasets)

      return(sim_output)

    }else if(advanced_dist == 'omic.two'){
      # Code to execute if advanced_dist is "omic.two"

      # features
      n_features_one <- vector_features[1]
      n_features_two <- vector_features[2]

      # Noise in the first detasets
      sigmas <- sigmas_vector[1]

      #num_factor = num_factors
      omic.one <- list()
      omic.two <- list()

      assigned_indices_samples <- divide_samples(n_samples = n_samples, num=n_factors, min_size=10)
      #print(assigned_indices_samples)

      # sample indices
      # Check the max number of factors in assigned_indices_samples
      max_factors <- length(assigned_indices_samples)# Same as n_factors

      # Generate random alpha values based on the max_factors
      list_alphas <- list()
      for (i in 1:max_factors) {
        list_alphas[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
      }

      # Assign corresponding values to alpha variables based on assigned_indices_samples
      for (i in seq_along(assigned_indices_samples)) {
        indices <- assigned_indices_samples[[i]]
        list_alphas[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
      }
      list_gammas <- list_alphas
      # Add the separated list of the indices selected for omics of multiple features
      # Initial vector
      vector <- c(1:max_factors)

      # Specify the number of elements to sample for the shared vector
      num_shared <- 0#ceiling(runif(1, min = 1, max = length(vector) - 1))

      # Select the 'shared' elements
      if (num_shared > 0) {
        shared <- sample(vector, num_shared)
        # Remove the selected shared elements from the vector
        remain_vector <- vector[!vector %in% shared]
      }else{
        shared <- 0
        remain_vector <- vector
      }

      # Shuffle the remaining vector
      if(length(remain_vector == 1)){
        shuffled_remain = remain_vector
      } else{
        shuffled_remain <- sample(remain_vector)
      }
      # Split the shuffled elements equally or nearly equally between omic_one_unique and omic_two_unique
      num_elements <- length(shuffled_remain)
      split_point <- sample(0:num_elements, 1) # Randomly select a split point

      omic_one_unique <- 0#shuffled_remain[0:split_point] # First part goes to omic_one_unique
      omic_two_unique <- remain_vector#shuffled_remain[!shuffled_remain %in% omic_one_unique] # Second part goes to omic_two_unique

      # Print the results
      list(
        shared = shared,
        omic_one_unique = omic_one_unique,
        omic_two_unique = omic_two_unique
      )

      # End of indices assignment

      # Assignment of factors
      list_omic_one_factors = c(omic_one_unique)
      list_omic_two_factors = c(omic_two_unique)
      #all_indices = c(shared, omic_one_unique, omic_two_unique)

      factor_xtics <- list(
        shared = shared,
        omic_one_unique = omic_one_unique,
        omic_two_unique = omic_two_unique
      )

      ### * * * * * * * * * * * * * * * * * * * * ** * * * * First OMIC data * * * * * * * * * * * * * * * * * * * * * * * * * *#

      # FeaturesD signal start_and_end
      assigned_indices_features = feature_selection_one(n_features_one = n_features_one, num.factor = "multiple", no_factor = n_factors)

      # Generate random beta values based on the max_factors
      list_betas <- list()

      # Create only one vector in list_betas using the random index
      for (i in seq_along(assigned_indices_features)) {
        list_betas[[i]] <- rnorm(n_features_one, 0, 0.05)
      }

      # Loop through the ordered indices in list_omic_one_factors
      for (i in seq_along(assigned_indices_features)) {
        indices_ns <- assigned_indices_features[[i]]  # Get the corresponding indices from assigned_indices_features
        if (length(indices_ns) > 0) {
          # Initialize list_betas[[i]] if it's not already
          if (is.null(list_betas[[i]])) {
            list_betas[[i]] <- numeric(max(indices_ns))
          }
          # Assign values to the appropriate indices in list_betas
          list_betas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = 0, sd = 0.05)
        }
      }

      # Renaming list_betas to desired factor number
      list_omic_factors <- seq_along(1:n_factors)
      for (i in seq_along(list_omic_factors)) {
        names(list_betas)[i] <- paste0("beta", list_omic_factors[[i]])#list_omic_one_factors[[i]]
      }

      pattern_alpha <-"^alpha\\d+"
      matches_alpha <- grepl(pattern_alpha, names(list_alphas), perl = TRUE)
      n_alpha = length(matches_alpha[matches_alpha == TRUE])

      pattern_beta <- "^beta\\d+"
      matches_beta <- grepl(pattern_beta, names(list_betas), perl = TRUE)
      n_beta = length(matches_beta[matches_beta == TRUE])

      # Initialize the data list to store the results
      data_list_i <- list()

      # Sort the list by names in ascending order
      list_betas <- list_betas[order(names(list_betas))]
      list_alphas <- list_alphas[order(names(list_alphas))]

      # Extract the last numeric values from the names of the vectors
      alphas_names <- as.numeric(gsub("[^0-9]", "", names(list_alphas)))
      betas_names <- as.numeric(gsub("[^0-9]", "", names(list_betas)))

      # Find common numbers between the last values of the vector names
      common_names <- intersect(alphas_names, betas_names)

      # Filter the vectors in each list based on the common names
      list_alphas <- list_alphas[paste0("alpha", common_names)]
      list_betas <- list_betas[paste0("beta", common_names)]

      # Loop through each alpha and beta combination
      for (i in 1:min(length(list_alphas), length(list_betas))) {
        data_i <- list_alphas[[i]] %*% t(list_betas[[i]][1:n_features_one])
        data_list_i[[paste0("data.", i)]] <- data_i
      }

      # Combine the results into a single data variable
      data.1 <- Reduce(`+`, data_list_i)

      eps1 <- rnorm(n_samples * n_features_one, 0, sigmas) # noise
      omic.one <- list()
      omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
      colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
      rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

      #dataset <- omic.one[[1]]
      #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")

      ### * * * * * * * * * * * * * * * * * * * * * * * Second OMIC data * * * * * * * * * * * * * * * * * * * * * *

      # Noise in the first detasets
      sigmas <- sigmas_vector[2]

      # Generate random gamma values based on the max_factors
      list_gammas <- list_alphas
      # Replace 'alpha' with 'gamma' in the names of list_gammas
      names(list_gammas) <- gsub("alpha", "gamma", names(list_gammas))

      # Extract the numeric parts from the names of list_gammas
      numeric_part <- as.numeric(gsub("\\D", "", names(list_gammas)))

      # Retain elements where the numeric part of the name matches values in list_omic_two_factors
      list_gammas <- list_gammas[numeric_part %in% list_omic_factors]

      # Features
      assigned_indices_features_omic.two = feature_selection_two(n_features_two = n_features_two, num.factor = "multiple", no_factor = length(list_omic_factors))
      (shuffled_assigned_indices_features_omic.two <- sample(assigned_indices_features_omic.two))
      # Empty list
      list_deltas <- list()

      # Create only one vector in list_deltas using the random index
      for (i in seq_along(shuffled_assigned_indices_features_omic.two)) {
        list_deltas[[i]] <- rnorm(n_features_two, 0, 0.05)
      }

      for (i in seq_along(shuffled_assigned_indices_features_omic.two)) {
        indices_ns <- shuffled_assigned_indices_features_omic.two[[i]]  # Get the corresponding indices from assigned_indices_features_omic.two
        # Remove NA values from indices_ns
        indices_ns <- indices_ns[!is.na(indices_ns)]
        if (length(indices_ns) > 0) {
          # Initialize list_deltas[[i]] if it's not already
          if (is.null(list_deltas[[i]])) {
            list_deltas[[i]] <- numeric(max(indices_ns, na.rm = TRUE))  # Use na.rm = TRUE to handle potential NA
          }
          # Ensure indices are within the range of list_deltas[[i]]
          if (max(indices_ns) > length(list_deltas[[i]])) {
            length(list_deltas[[i]]) <- max(indices_ns)
          }
          # Assign values to the appropriate indices in list_deltas
          list_deltas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (5.0 + 0.4 * i), sd = 0.05)
        }
      }

      # Renaming list_deltas to desired factor number
      list_omic_factors <- seq_along(1:n_factors)
      for (i in seq_along(list_omic_factors)) {
        names(list_deltas)[i] <- paste0("delta", list_omic_factors[[i]])#list_omic_one_factors[[i]]
      }

      pattern_gamma <-"^gamma\\d+"
      matches_gamma <- grepl(pattern_gamma, names(list_gammas), perl = TRUE)
      n_gamma = length(matches_alpha[matches_gamma == TRUE])

      pattern_delta <- "^delta\\d+"
      matches_delta <- grepl(pattern_delta, names(list_deltas), perl = TRUE)
      n_delta = length(matches_delta[matches_delta == TRUE])

      # Initialize the data list to store the results
      data_list_j <- list()

      # Sort the list by names in ascending order
      list_gammas <- list_gammas[order(names(list_gammas))]
      list_deltas <- list_deltas[order(names(list_deltas))]

      # Extract the last numeric values from the names of the vectors
      gammas_names <- as.numeric(gsub("[^0-9]", "", names(list_gammas)))
      deltas_names <- as.numeric(gsub("[^0-9]", "", names(list_deltas)))

      # Find common numbers between the last values of the vector names
      common_names <- intersect(gammas_names, deltas_names)

      # Filter the vectors in each list based on the common names
      list_gammas <- list_gammas[paste0("gamma", common_names)]
      list_deltas <- list_deltas[paste0("delta", common_names)]

      # Loop through each alpha and beta combination
      for (j in 1:min(length(list_gammas), length(list_deltas))) {
        data_j <- list_gammas[[j]] %*% t(list_deltas[[j]][1:n_features_two])
        data_list_j[[paste0("data.", j)]] <- data_j
      }

      # Combine the results into a single data variable
      data.2 <- Reduce(`+`, data_list_j)

      eps2 <- rnorm(n_samples * n_features_two, 0, sigmas)
      omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
      colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
      rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

      #dataset <- data.2#dataset <- omic.two[[1]]
      #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")

      simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

      concatenated_datasets <- list()
      for (i in 1:length(simulated_datasets$object.two)) {
        concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
        concatenated_datasets[[i]] <- concatenated_data
      }

      # sim_output <- list(concatenated_datasets = concatenated_datasets, list_alphas = list_alphas, list_betas = list_betas, list_gammas = list_gammas, list_deltas = list_deltas)#simulated_datasets=simulated_datasets)
      sim_output <- list(concatenated_datasets=concatenated_datasets, omic.one=omic.one, omic.two=omic.two, assigned_indices_samples = assigned_indices_samples, assigned_indices_features = assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, factor_xtics=factor_xtics, list_alphas=list_alphas, list_gammas=list_gammas, list_betas=list_betas, list_deltas=list_deltas)#simulated_datasets=simulated_datasets)

      return(sim_output)

    }else if(advanced_dist == 'exclusive'){
      # Code to execute if advanced_dist is "exclusive"

      # features
      n_features_one <- vector_features[1]
      n_features_two <- vector_features[2]

      # Noise in the first detasets
      sigmas <- sigmas_vector[1]

      #num_factor = num_factors
      omic.one <- list()
      omic.two <- list()

      assigned_indices_samples <- divide_samples(n_samples = n_samples, num=n_factors, min_size=10)
      #print(assigned_indices_samples)

      # sample indices
      # Check the max number of factors in assigned_indices_samples
      max_factors <- length(assigned_indices_samples)# Same as n_factors

      # Generate random alpha values based on the max_factors
      list_alphas <- list()
      for (i in 1:max_factors) {
        list_alphas[[paste0("alpha", i)]] <- rnorm(n_samples, 0, 0.05)
      }

      # Assign corresponding values to alpha variables based on assigned_indices_samples
      for (i in seq_along(assigned_indices_samples)) {
        indices <- assigned_indices_samples[[i]]
        list_alphas[[i]][indices] <- rnorm(length(indices), (3 + 0.5*i), 0.05)  # Adjust values dynamically
      }
      list_gammas <- list_alphas
      # Add the separated list of the indices selected for omics of multiple features
      # Initial vector
      vector <- c(1:max_factors)

      # Specify the number of elements to sample for the shared vector
      num_shared <- 0#ceiling(runif(1, min = 1, max = length(vector) - 1))

      # Select the 'shared' elements
      if (num_shared > 0) {
        shared <- sample(vector, num_shared)
        # Remove the selected shared elements from the vector
        remain_vector <- vector[!vector %in% shared]
      }else{
        shared <- 0
        remain_vector <- vector
      }

      # Shuffle the remaining vector
      if(length(remain_vector == 1)){
        shuffled_remain = remain_vector
      } else{
        shuffled_remain <- sample(remain_vector)
      }
      # Split the shuffled elements equally or nearly equally between omic_one_unique and omic_two_unique
      num_elements <- length(shuffled_remain)
      split_point <- sample(1:(num_elements-1), 1) # Randomly select a split point

      omic_one_unique <- shuffled_remain[1:split_point] # First part goes to omic_one_unique
      omic_two_unique <- shuffled_remain[!shuffled_remain %in% omic_one_unique] # Second part goes to omic_two_unique

      # Print the results
      list(
        shared = shared,
        omic_one_unique = omic_one_unique,
        omic_two_unique = omic_two_unique
      )

      # End of indices assignment

      # Assignment of factors
      list_omic_one_factors = omic_one_unique
      list_omic_two_factors = omic_two_unique
      #all_indices = c(shared, omic_one_unique, omic_two_unique)

      factor_xtics <- list(
        shared = shared,
        omic_one_unique = omic_one_unique,
        omic_two_unique = omic_two_unique
      )

      ### * * * * * * * * * * * * * * * * * * * * ** * * * * First OMIC data * * * * * * * * * * * * * * * * * * * * * * * * * *#

      # FeaturesD signal start_and_end
      assigned_indices_features = feature_selection_one(n_features_one = n_features_one, num.factor = "multiple", no_factor = length(list_omic_one_factors))

      # Generate random beta values based on the max_factors
      list_betas <- list()

      # Shuffle the list
      (shuffled_assigned_indices_features <- assigned_indices_features)# sample(assigned_indices_features)) # Do not reshuffle first

      # Create only one vector in list_betas using the random index
      for (i in seq_along(assigned_indices_features)) {
        list_betas[[i]] <- rnorm(n_features_one, 0, 0.05)
      }

      # Loop through the ordered indices in list_omic_one_factors
      for (i in seq_along(assigned_indices_features)) {
        indices_ns <- assigned_indices_features[[i]]  # Get the corresponding indices from assigned_indices_features
        if (length(indices_ns) > 0) {
          # Initialize list_betas[[i]] if it's not already
          if (is.null(list_betas[[i]])) {
            list_betas[[i]] <- numeric(max(indices_ns))
          }
          # Assign values to the appropriate indices in list_betas
          list_betas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (5.0 + 0.4 * i), sd = 0.05)
        }
      }

      # Renaming list_betas to desired factor number
      for (i in seq_along(list_omic_one_factors)) {
        names(list_betas)[i] <- paste0("beta", list_omic_one_factors[[i]])#list_omic_one_factors[[i]]
      }

      pattern_alpha <-"^alpha\\d+"
      matches_alpha <- grepl(pattern_alpha, names(list_alphas), perl = TRUE)
      n_alpha = length(matches_alpha[matches_alpha == TRUE])

      pattern_beta <- "^beta\\d+"
      matches_beta <- grepl(pattern_beta, names(list_betas), perl = TRUE)
      n_beta = length(matches_beta[matches_beta == TRUE])

      # Initialize the data list to store the results
      data_list_i <- list()

      # Sort the list by names in ascending order
      list_betas <- list_betas[order(names(list_betas))]
      list_alphas <- list_alphas[order(names(list_alphas))]

      # Extract the last numeric values from the names of the vectors
      alphas_names <- as.numeric(gsub("[^0-9]", "", names(list_alphas)))
      betas_names <- as.numeric(gsub("[^0-9]", "", names(list_betas)))

      # Find common numbers between the last values of the vector names
      common_names <- intersect(alphas_names, betas_names)

      # Filter the vectors in each list based on the common names
      list_alphas <- list_alphas[paste0("alpha", common_names)]
      list_betas <- list_betas[paste0("beta", common_names)]

      # Loop through each alpha and beta combination
      for (i in 1:min(length(list_alphas), length(list_betas))) {
        data_i <- list_alphas[[i]] %*% t(list_betas[[i]][1:n_features_one])
        data_list_i[[paste0("data.", i)]] <- data_i
      }

      # Combine the results into a single data variable
      data.1 <- Reduce(`+`, data_list_i)

      eps1 <- rnorm(n_samples * n_features_one, 0, sigmas) # noise
      omic.one <- list()
      omic.one[[length(omic.one) + 1]] <- matrix(data.1, n_samples, n_features_one) + matrix(eps1, n_samples, n_features_one) # signal + noise
      colnames(omic.one[[length(omic.one)]]) <- paste0('omic1_feature_', seq_len(n_features_one))
      rownames(omic.one[[length(omic.one)]]) <- paste0('sample_', seq_len(n_samples))

      #dataset <- omic.one[[1]]
      #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")

      ### * * * * * * * * * * * * * * * * * * * * ** * * * * Second OMIC data * * * * * * * * * * * * * * * * * * * * * * * * * *

      # Noise in the first detasets
      sigmas <- sigmas_vector[2]

      # Generate random gamma values based on the max_factors
      #list_gammas <- list_alphas
      # Replace 'alpha' with 'gamma' in the names of list_gammas
      names(list_gammas) <- gsub("alpha", "gamma", names(list_gammas))

      # Extract the numeric parts from the names of list_gammas
      numeric_part <- as.numeric(gsub("\\D", "", names(list_gammas)))

      # Retain elements where the numeric part of the name matches values in list_omic_two_factors
      list_gammas <- list_gammas[numeric_part %in% list_omic_two_factors]

      # Features
      assigned_indices_features_omic.two = feature_selection_two(n_features_two = n_features_two, num.factor = "multiple", no_factor = length(list_omic_two_factors))

      # Empty list
      list_deltas <- list()

      # Create only one vector in list_deltas using the random index
      for (i in seq_along(assigned_indices_features_omic.two)) {
        list_deltas[[i]] <- rnorm(n_features_two, 0, 0.05)
      }

      # # Loop through the ordered indices in list_omic_one_factors
      for (i in seq_along(assigned_indices_features_omic.two)) {
        indices_ns <- assigned_indices_features_omic.two[[i]]  # Get the corresponding indices from assigned_indices_features_omic.two
        # Remove NA values from indices_ns
        indices_ns <- indices_ns[!is.na(indices_ns)]
        if (length(indices_ns) > 0) {
          # Initialize list_deltas[[i]] if it's not already
          if (is.null(list_deltas[[i]])) {
            list_deltas[[i]] <- numeric(max(indices_ns, na.rm = TRUE))  # Use na.rm = TRUE to handle potential NA
          }
          # Ensure indices are within the range of list_deltas[[i]]
          if (max(indices_ns) > length(list_deltas[[i]])) {
            length(list_deltas[[i]]) <- max(indices_ns)
          }
          # Assign values to the appropriate indices in list_deltas
          list_deltas[[i]][indices_ns] <- rnorm(length(indices_ns), mean = (5.0 + 0.4 * i), sd = 0.05)
        }
      }

      # Renaming list_deltas to desired factor number
      for (i in seq_along(list_omic_two_factors)) {
        names(list_deltas)[i] <- paste0("delta", list_omic_two_factors[[i]])#list_omic_one_factors[[i]]
      }

      pattern_gamma <-"^gamma\\d+"
      matches_gamma <- grepl(pattern_gamma, names(list_gammas), perl = TRUE)
      n_gamma = length(matches_alpha[matches_gamma == TRUE])

      pattern_delta <- "^delta\\d+"
      matches_delta <- grepl(pattern_delta, names(list_deltas), perl = TRUE)
      n_delta = length(matches_delta[matches_delta == TRUE])

      # Initialize the data list to store the results
      data_list_j <- list()

      # Sort the list by names in ascending order
      list_gammas <- list_gammas[order(names(list_gammas))]
      list_deltas <- list_deltas[order(names(list_deltas))]

      # Extract the last numeric values from the names of the vectors
      gammas_names <- as.numeric(gsub("[^0-9]", "", names(list_gammas)))
      deltas_names <- as.numeric(gsub("[^0-9]", "", names(list_deltas)))

      # Find common numbers between the last values of the vector names
      common_names <- intersect(gammas_names, deltas_names)

      # Filter the vectors in each list based on the common names
      list_gammas <- list_gammas[paste0("gamma", common_names)]
      list_deltas <- list_deltas[paste0("delta", common_names)]

      # Loop through each alpha and beta combination
      for (j in 1:min(length(list_gammas), length(list_deltas))) {
        data_j <- list_gammas[[j]] %*% t(list_deltas[[j]][1:n_features_two])
        data_list_j[[paste0("data.", j)]] <- data_j
      }

      # Combine the results into a single data variable
      data.2 <- Reduce(`+`, data_list_j)

      eps2 <- rnorm(n_samples * n_features_two, 0,sigmas)
      omic.two[[length(omic.two) + 1]] <- matrix(data.2, n_samples, n_features_two) + matrix(eps2, n_samples, n_features_two) # signal + noise
      colnames(omic.two[[length(omic.two)]]) <- paste0('omic2_feature_', seq_len(n_features_two))
      rownames(omic.two[[length(omic.two)]]) <- paste0('sample_', seq_len(n_samples))

      #dataset <- data.2#dataset <- omic.two[[1]]
      #image(c(1:dim(dataset)[1]), c(1:dim(dataset)[2]), dataset, ylab = "Features", xlab = "Samples")

      simulated_datasets <- list(object.one = omic.one, object.two = omic.two)

      concatenated_datasets <- list()
      for (i in 1:length(simulated_datasets$object.two)) {
        concatenated_data <- cbind(simulated_datasets$object.two[[i]],simulated_datasets$object.one[[i]])
        concatenated_datasets[[i]] <- concatenated_data
      }

      # sim_output <- list(concatenated_datasets = concatenated_datasets, list_alphas = list_alphas, list_betas = list_betas, list_gammas = list_gammas, list_deltas = list_deltas)#simulated_datasets=simulated_datasets)
      sim_output <- list(concatenated_datasets=concatenated_datasets, omic.one=omic.one, omic.two=omic.two, assigned_indices_samples = assigned_indices_samples, assigned_indices_features = assigned_indices_features, assigned_indices_features_omic.two=assigned_indices_features_omic.two, factor_xtics=factor_xtics, list_alphas=list_alphas, list_gammas=list_gammas, list_betas=list_betas, list_deltas=list_deltas)#simulated_datasets=simulated_datasets)

      return(sim_output)
    }
  }else{
      message <- "Error: If 'num.factor' must be provided with 'single' or 'multiple'"
      print(message) # Display message immediately
      return(message)
  }
}

# Testing
#set.seed(123)
#output_sim <- OmixCraftHD(vector_features = c(2000,3000), sigmas_vector=c(8,5), n_samples=100,
#                          n_factors=5, num.factor='multiple', advanced_dist='mixed')
#output_sim$factor_xtics
#plot_simData(output_sim, type = 'heatmap')
#plot_factor(1)
