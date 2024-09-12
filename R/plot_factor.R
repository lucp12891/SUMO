#' @name plot_factor
#' @title Visualization of factor scores
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point after_stat element_blank labs scale_color_viridis_c theme_minimal  theme_bw ggtitle theme element_text scale_fill_viridis_c
#' @param sim_object R object containing data to be plotted
#' @param factor_num Factor to be plotted.
#' @export
  plot_factor <- function(sim_object = NULL, factor_num = NULL) {
    # libraries
    #library(ggplot2); require(ggplot2)
    #library(gridExtra); require(gridExtra)

    #Run if no parameter is passed
    if (is.null(sim_object) || is.null(factor_num)) {
      vector_features = c(2000,2000)
      n_samples = 50
      sigmas_vector = c(3,5)
      n_factors = 3
      num.factor = 'multiple'
      advanced_dist = NULL

      sim_object <- OmixCraftHD(vector_features = c(2000,2000), n_samples = 50, sigmas_vector = c(3,5), n_factors = 3, num.factor = 'multiple', advanced_dist = NULL)
      factor_num = 1

      # Provide multi-line feedback to the user about the generated data and parameters
      message("Note:")
      message(" - No parameters were passed to the function.")
      message(" - Simulation using the default settings.")
      message(paste(" - vector_features:", vector_features))
      message(paste(" - n_samples:", n_samples))
      message(paste(" - sigmas_vector:", sigmas_vector))
      message(paste(" - n_factors:", n_factors))
      message(paste(" - num.factor:", num.factor))
      message(paste(" - advanced_dist:", advanced_dist))
      message(paste(" - Default factor_num:", factor_num))

    }else if(!is.null(sim_object) && is.null(factor_num)){
      sim_object = sim_object
      factor_num = 1
      message("Note:")
      message(" - Number of factor not passed to the function.")
      message(paste(" - Default factor_num:", sim_object))
      message(paste(" - Default factor_num:", factor_num))
    }
    # Load necessary libraries
    #library(ggplot2::ggplot2)
    #library(gridExtra)

    # Combine the two lists
    combined_list <- c(sim_object$list_alphas, sim_object$list_gammas)

    # Extract the names of the combined list
    combined_names <- names(combined_list)

    # Extract the digits at the end of each name using a regular expression
    # This will capture the numeric part at the end of each name
    end_digits <- sub(".*?(\\d+)$", "\\1", combined_names)

    # Create a named vector to track unique elements by their ending digit
    unique_end_digits <- !duplicated(end_digits)

    # Create the final list with only one occurrence for each digit
    factor_scores <- combined_list[unique_end_digits]

    # Replace 'alpha' and 'gamma' with score in the names of factor_scores
    names(factor_scores) <- gsub("alpha", "score", names(factor_scores))
    names(factor_scores) <- gsub("gamma", "score", names(factor_scores))

    # Extract data from the sim_object
    #data <- (output_sim$concatenated_datasets[[1]])
    data <- data.frame(sim_object$concatenated_datasets[[1]])
    samples <- rownames(data)

    # Check if num is in the range of theta indices
    if (factor_num %in% seq_along(factor_scores)) {

      factor <- factor_scores[[paste0('score', factor_num)]] # Retrieve the corresponding theta value

      factor_df <- data.frame(samples = samples, factor = factor)       # Create a data frame with names and retrieved values

      # Extract numeric values using regular expressions
      factor_df$Index <- as.numeric(gsub("\\D", "", factor_df$samples))

      # Sort the data frame based on Iteration and Sigma
      factor_df_sorted <- factor_df[order(factor_df$Index), ]

      # scatter plot
      if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x=Index, y=factor, color = factor)) +
        geom_point()+
        labs(x = "Samples", y = "Factor score", color = paste("Magnitude")) +
        #theme(legend.position="none")+
        scale_color_viridis_c("Magnitude") +
        theme_bw()
      } else {
        stop("ggplot2 package is required but not installed.")
      }
      # return
      return(plot)

    } else if (factor_num == 'all'){
      #sim_object = output_sim
      data <- data.frame(sim_object$concatenated_datasets[[1]])
      samples <- rownames(data)

      factors <- names(c(factor_scores[grep("score", names(factor_scores))]))

      # Initialize plot_list to store plots
      plot_list <- list()

      # Create plots for each factor and store in the plot_list
      for (i in seq_along(factors)) {
        factor_name <- paste("score", i, sep = "")
        factor <- factor_scores[[factor_name]]
        factor_df <- data.frame(samples = samples, factor = factor)

        # Extract numeric values using regular expressions
        factor_df$Index <- as.numeric(gsub("\\D", "", factor_df$samples))

        # Sort the data frame based on Index
        factor_df_sorted <- factor_df[order(factor_df$Index), ]

        # Scatter plot
        if (requireNamespace("ggplot2", quietly = TRUE)) {
        plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x = Index, y = factor, color = factor)) +
          geom_point() +
          labs(x = "Samples", y = paste("Factor", i, "score"), color = "") +
          scale_color_viridis_c("") + #scale_color_viridis_c("Magnitude")
          theme_bw()
        } else {
          stop("ggplot2 package is required but not installed.")
        }
        # Store the plot in plot_list
        plot_list[[i]] <- plot
      }
      if (requireNamespace("gridExtra", quietly = TRUE)) {
      plot_arrange <- gridExtra::grid.arrange(grobs = plot_list, ncol = length(factors))
      } else {
        stop("The 'gridExtra' package is required but not installed. Please install it to use this functionality.")
      }
    } else {
      print("Provide the correct input.")
      return(NULL)
    }
  }
  # Suppressing the global variable warning
  utils::globalVariables(c("Index"))

  #################################################################################

# Examples
## plot_factor(sim_object, factor_num = 'all')
## plot_factor(output_sim, factor_num = 1)
## plot_factor(output_sim, factor_num = 5)
## plot_factor()
## plot_factor(output_sim)
