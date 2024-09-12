#' @name plot_weights
#' @title Visualizing the loading of the features
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point after_stat element_blank labs scale_color_viridis_c theme_minimal  theme_bw ggtitle theme element_text scale_fill_viridis_c
#' @importFrom gridExtra grid.arrange
#' @param sim_object R object containing data to be plotted
#' @param factor_num Factor to be plotted.
#' @param data Section of the integrated data to be plotted, omic.one or omic.two are the options
#' @param type Type of plot. Scatter plot and histogram are the only allowed plots
#' @importFrom rlang .data
#' @export
plot_weights <- function(sim_object = NULL, factor_num = 1, data = 'omic.one', type = 'scatter') {
  # Load necessary libraries
  #library(ggplot2); require(ggplot2)
  #library(gridExtra); require(gridExtra)

  #Run if no parameter is passed
  if (is.null(sim_object) || is.null(factor_num) ||is.null(data) ||is.null(type)) {
    vector_features = c(2000,2000)
    n_samples = 50
    sigmas_vector = c(3,5)
    n_factors = 3
    num.factor = 'multiple'
    type = 'scatter'
    data = 'omic.one'
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
    message(paste(" - Default data:", data))
    message(paste(" - Default type:", type))

  }

  if(!is.null(sim_object) && is.null(factor_num)){
    sim_object = sim_object
    factor_num = 1
    type = 'scatter'
    data = 'omic.one'
    message("Note:")
    message(" - Number of factor not passed to the function.")
    message(paste(" - Default factor_num:", sim_object))
    message(paste(" - Default factor_num:", factor_num))
    message(paste(" - Default data:", data))
    message(paste(" - Default type:", type))
  }


  if(type == 'scatter'){
  # Check if num is in the range of theta indices
  if(data == 'omic.one'){#}

    sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

    # Use grep to find column names starting with "omic1_feature_"
    features <- grep("^omic1_feature_", colnames(sim_data), value = TRUE)

    # Extract the last numeric values from the names of the vectors
    betas_names <- as.numeric(gsub("[^0-9]", "", names(sim_object$list_betas)))

    if (factor_num %in% betas_names) {

        factor <- sim_object$list_betas[[paste0('beta', factor_num)]] # Retrieve the corresponding theta value

        factor_df <- data.frame(features = features, factor = factor)       # Create a data frame with names and retrieved values

        # Extract numeric values using regular expressions
        factor_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_df$features))

        # Sort the data frame based on Iteration and Sigma
        factor_df_sorted <- factor_df[order(factor_df$Index), ]

        # scatter plot
        if (requireNamespace("ggplot2", quietly = TRUE)) {
        plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x=Index, y=factor, color = factor)) +
          geom_point()+
          labs(x = "Features", y = "Weights", color = paste("")) + #color = paste("Magnitude"))
          #theme(legend.position="none")+
          scale_color_viridis_c("") + #scale_color_viridis_c("Magnitude")
          theme_bw() +
          ggtitle("My Plot Title")
        } else {
          stop("ggplot2 package is required but not installed.")
        }
        # return
        return(plot)

      } else if (factor_num == 'all'){

        #data <- data.frame(sim_object$concatenated_datasets[[1]])
        #features <- colnames(data)

        factors <- names(c(sim_object$list_betas[grep("beta", names(sim_object$list_betas))]))

        # Initialize plot_list to store plots
        plot_list <- list()

        # Create plots for each factor and store in the plot_list
        for (i in seq_along(factors)) {
          factor_name <- paste("beta", i, sep = "")
          factor <- sim_object$list_betas[[factor_name]]
          factor_df <- data.frame(features = features, factor = factor)

          # Extract numeric values using regular expressions
          factor_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_df$features))

          # Sort the data frame based on Index
          factor_df_sorted <- factor_df[order(factor_df$Index), ]

          # Scatter plot
          if (requireNamespace("ggplot2", quietly = TRUE)) {
          plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x = Index, y = factor, color = factor)) +
            geom_point() +
            labs(x = "Features", y = paste("Factor_", i, "_Weights"), color = "") + #color = paste("Magnitude"))
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
        message <- paste0("Provide the correct input. Available options: ", paste(betas_names, collapse = ", "), " and all") # Message
        print(message) # Print the message and number of available elements
        return(message) # Return the message
      }
  } else if(data == 'omic.two'){
    sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

    # Use grep to find column names starting with "omic2_feature_"
    features <- grep("^omic2_feature_", colnames(sim_data), value = TRUE)

    # Extract the last numeric values from the names of the vectors
    deltas_names <- as.numeric(gsub("[^0-9]", "", names(sim_object$list_deltas)))

    if (factor_num %in% deltas_names) {
        factor <- sim_object$list_deltas[[paste0('delta', factor_num)]] # Retrieve the corresponding theta value

        factor_df <- data.frame(features = features, factor = factor)       # Create a data frame with names and retrieved values

        # Extract numeric values using regular expressions
        factor_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_df$features))

        # Sort the data frame based on Iteration and Sigma
        factor_df_sorted <- factor_df[order(factor_df$Index), ]

        # scatter plot
        if (requireNamespace("ggplot2", quietly = TRUE)) {
        plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x=Index, y=factor, color = factor)) +
          geom_point()+
          labs(x = "Features", y = "Weights", color = paste("")) +  #color = paste("Magnitude"))
          #theme(legend.position="none")+
          scale_color_viridis_c("") + #scale_color_viridis_c("Magnitude")
          theme_bw()
        } else {
          stop("ggplot2 package is required but not installed.")
        }
        # return
        return(plot)

    } else if (factor_num == 'all'){
      #data <- data.frame(sim_object$concatenated_datasets[[1]])
      #features <- colnames(data)

      factors <- deltas_names #names(c(sim_object$list_deltas[grep("delta", names(sim_object$list_deltas))]))

      # Initialize plot_list to store plots
      plot_list <- list()

      # Create plots for each factor and store in the plot_list
      for (i in factors) {
        factor_name <- paste("delta", i, sep = "")
        factor <- sim_object$list_deltas[[factor_name]]
        factor_df <- data.frame(features = features, factor = factor)

        # Extract numeric values using regular expressions
        factor_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_df$features))

        # Sort the data frame based on Index
        factor_df_sorted <- factor_df[order(factor_df$Index), ]

        # Scatter plot
        if (requireNamespace("ggplot2", quietly = TRUE)) {
        plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x = Index, y = factor, color = factor)) +
          geom_point() +
          labs(x = "Features", y = paste("Factor", i, "Weights"), color = "") + #color = "Magnitude")
          scale_color_viridis_c("") + scale_color_viridis_c("Magnitude")
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
      message <- paste0("Provide the correct input. Available options: ", paste(deltas_names, collapse = ", "), " and all") # Message
      #print(message) # Print the message and number of available elements
      return(message) # Return the message
    }
  } else if(data == 'integrated'){
    #  # list_betas # list_deltas
    #  # merge the lists
    data <- data.frame(sim_object$concatenated_datasets[[1]])
    features <- colnames(data)

    # Initialize empty vectors to store beta names and numeric values
    list_betas <- sim_object$list_betas
    list_deltas <- sim_object$list_deltas

    # Combine vectors as rows into a data frame
    betas_df <- data.frame(do.call(cbind, list_betas))
    deltas_df <- data.frame(do.call(cbind, list_deltas))

    num_cols_beta <- ncol(betas_df)
    num_cols_delta <- ncol(deltas_df)

    # Extract the last numeric values from the names of the vectors
    betas_names <- as.numeric(gsub("[^0-9]", "", names(sim_object$list_betas)))
    deltas_names <- as.numeric(gsub("[^0-9]", "", names(sim_object$list_deltas)))

    #new_col_beta_names <- paste0("factor", seq_len(num_cols_beta))
    #new_col_delta_names <- paste0("factor", seq_len(num_cols_delta))
    new_col_beta_names <- paste0("factor", betas_names)
    new_col_delta_names <- paste0("factor", deltas_names)

    # Assign new column names
    colnames(betas_df) <- new_col_beta_names
    colnames(deltas_df) <- new_col_delta_names

    # Omic1
    if (factor_num %in% betas_names && factor_num %in% deltas_names) {#if (factor_num %in% seq_along(sim_object$list_betas)) {

      sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

      # Use grep to find column names starting with "omic1_feature_"
      features.omic1 <- grep("^omic1_feature_", colnames(sim_data), value = TRUE)

      factor_omic1 <- sim_object$list_betas[[paste0('beta', factor_num)]] # Retrieve the corresponding theta value

      factor_omic1_df <- data.frame(features = features.omic1, factor = factor_omic1)       # Create a data frame with names and retrieved values

      # Extract numeric values using regular expressions
      factor_omic1_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_omic1_df$features))

      # Sort the data frame based on Iteration and Sigma
      factor_omic1_df_sorted <- factor_omic1_df[order(factor_omic1_df$Index), ]

      # scatter plot
      if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot.a <- ggplot2::ggplot(data = factor_omic1_df_sorted, aes(x=Index, y=factor, color = factor)) +
        geom_point()+
        labs(x = "Features", y = "Weights", color = paste("")) + #color = "Magnitude")
        theme(
          plot.title = element_text(size = 8, hjust = 0.5)  # Adjust title size and center it
        ) +
        scale_color_viridis_c("") + #scale_color_viridis_c("Magnitude")
        theme_bw()+
        ggtitle(paste("Data: Omic.one, Factor:", factor_num))
      } else {
        stop("ggplot2 package is required but not installed.")
      }

      sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

      # Use grep to find column names starting with "omic1_feature_"
      features_omic2 <- grep("^omic2_feature_", colnames(sim_data), value = TRUE)

      factor_omic2 <- sim_object$list_deltas[[paste0('delta', factor_num)]] # Retrieve the corresponding theta value

      factor_omic2_df <- data.frame(features = features_omic2, factor = factor_omic2)       # Create a data frame with names and retrieved values

      # Extract numeric values using regular expressions
      factor_omic2_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_omic2_df$features))

      # Sort the data frame based on Iteration and Sigma
      factor_omic2_df_sorted <- factor_omic2_df[order(factor_omic2_df$Index), ]

      # scatter plot
      if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot.b <- ggplot2::ggplot(data = factor_omic2_df_sorted, aes(x=Index, y=factor, color = factor)) +
        geom_point()+
        labs(x = "Features", y = "Weights", color = paste("")) + #color = paste("Magnitude"))
        theme(
          plot.title = element_text(size = 8, hjust = 0.5)  # Adjust title size and center it
        ) +
        scale_color_viridis_c("") + #scale_color_viridis_c("Magnitude")
        theme_bw()+
        ggtitle(paste("Data: Omic.two, Factor:", factor_num))
      } else {
        stop("ggplot2 package is required but not installed.")
      }
      # return
      plot_list <- list(plot.a, plot.b)
      if (requireNamespace("gridExtra", quietly = TRUE)) {
      plot_arrange <- gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
      } else {
        stop("The 'gridExtra' package is required but not installed. Please install it to use this functionality.")
      }
    }else if(factor_num %in% betas_names && !factor_num %in% deltas_names){
      # Omic1
        sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

        # Use grep to find column names starting with "omic1_feature_"
        features.omic1 <- grep("^omic1_feature_", colnames(sim_data), value = TRUE)

        factor_omic1 <- sim_object$list_betas[[paste0('beta', factor_num)]] # Retrieve the corresponding theta value

        factor_omic1_df <- data.frame(features = features.omic1, factor = factor_omic1)       # Create a data frame with names and retrieved values

        # Extract numeric values using regular expressions
        factor_omic1_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_omic1_df$features))

        # Sort the data frame based on Iteration and Sigma
        factor_omic1_df_sorted <- factor_omic1_df[order(factor_omic1_df$Index), ]

        # scatter plot
        if (requireNamespace("ggplot2", quietly = TRUE)) {
        plot.a <- ggplot2::ggplot(data = factor_omic1_df_sorted, aes(x=Index, y=factor, color = factor)) +
          geom_point()+
          labs(x = "Features", y = "Weights", color = paste("")) + #color = paste("Magnitude"))
          theme(
            plot.title = element_text(size = 8, hjust = 0.5)  # Adjust title size and center it
          ) +
          scale_color_viridis_c("") + #scale_color_viridis_c("Magnitude")
          theme_bw()+
          ggtitle(paste("Data: Omic.one, Factor:", factor_num))
        } else {
          stop("ggplot2 package is required but not installed.")
        }
        # return
        return(plot.a)

    }else if(!factor_num %in% betas_names && factor_num %in% deltas_names){
      # Omic2
      sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

      # Use grep to find column names starting with "omic1_feature_"
      features_omic2 <- grep("^omic2_feature_", colnames(sim_data), value = TRUE)

      factor_omic2 <- sim_object$list_deltas[[paste0('delta', factor_num)]] # Retrieve the corresponding theta value

      factor_omic2_df <- data.frame(features = features_omic2, factor = factor_omic2)       # Create a data frame with names and retrieved values

      # Extract numeric values using regular expressions
      factor_omic2_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_omic2_df$features))

      # Sort the data frame based on Iteration and Sigma
      factor_omic2_df_sorted <- factor_omic2_df[order(factor_omic2_df$Index), ]

      # scatter plot
      if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot.b <- ggplot2::ggplot(data = factor_omic2_df_sorted, aes(x=Index, y=factor, color = factor)) +
        geom_point()+
        labs(x = "Features", y = "Weights", color = paste("")) + #color = paste("Magnitude"))
        theme(
          plot.title = element_text(size = 8, hjust = 0.5)  # Adjust title size and center it
        ) +
        scale_color_viridis_c("") + #scale_color_viridis_c("Magnitude")
        theme_bw()+
        ggtitle(paste("Data: Omic.two, Factor:", factor_num))
      } else {
        stop("ggplot2 package is required but not installed.")
      }
      # return
      return(plot.b)

    }
  }
  } else if (type == 'histogram'){
   # Check if num is in the range of theta indices
   if(data == 'omic.one'){#}

     sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

     # Use grep to find column names starting with "omic1_feature_"
     features <- grep("^omic1_feature_", colnames(sim_data), value = TRUE)

     if (factor_num %in% betas_names) {#if (factor_num %in% seq_along(sim_object$list_betas)) {

       factor <- sim_object$list_betas[[paste0('beta', factor_num)]] # Retrieve the corresponding theta value

       factor_df <- data.frame(features = features, factor = factor)       # Create a data frame with names and retrieved values

       # Extract numeric values using regular expressions
       factor_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_df$features))

       # Sort the data frame based on Iteration and Sigma
       factor_df_sorted <- factor_df[order(factor_df$Index), ]

       # histogram plot
       if (requireNamespace("ggplot2", quietly = TRUE)) {
       plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x = factor, fill = after_stat(x))) +
         geom_histogram(bins=100) +  # Use geom_bar() instead of geom_histogram()
         labs(x = "Loadings", y = "Frequency", fill = "") + #fill = "Magnitude")
         scale_fill_viridis_c() +  # Use scale_fill_viridis_c() for fill color
         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         theme_minimal()
       } else {
         stop("ggplot2 package is required but not installed.")
       }
       # return
       return(plot)

     } else if (factor_num == 'all'){

       #data <- data.frame(sim_object$concatenated_datasets[[1]])
       #features <- colnames(data)

       factors <- betas_names#names(c(sim_object$list_betas[grep("beta", names(sim_object$list_betas))]))

       # Initialize plot_list to store plots
       plot_list <- list()

       # Create plots for each factor and store in the plot_list
       for (i in seq_along(factors)) {
         factor_name <- paste("beta", i, sep = "")
         factor <- sim_object$list_betas[[factor_name]]
         factor_df <- data.frame(features = features, factor = factor)

         # Extract numeric values using regular expressions
         factor_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_df$features))

         # Sort the data frame based on Index
         factor_df_sorted <- factor_df[order(factor_df$Index), ]

         # Scatter plot
         if (requireNamespace("ggplot2", quietly = TRUE)) {
         plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x = factor, fill = after_stat(x))) +
           geom_histogram(bins=100) +  # Use geom_bar() instead of geom_histogram()
           labs(x = "Loadings", y = "Frequency", fill = "") + #fill = "Magnitude")
           scale_fill_viridis_c() +  # Use scale_fill_viridis_c() for fill color
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
           theme_minimal()
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
       message <- paste0("Provide the correct input. Available options: ", paste(betas_names, collapse = ", "), " and all") # Message
       #print(message) # Print the message and number of available elements
       return(message) # Return the message

     }
   } else if(data == 'omic.two'){
     sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

     # Use grep to find column names starting with "omic2_feature_"
     features <- grep("^omic2_feature_", colnames(sim_data), value = TRUE)

     if (factor_num %in% deltas_names){ #if (factor_num %in% seq_along(sim_object$list_deltas)) {
       factor <- sim_object$list_deltas[[paste0('delta', factor_num)]] # Retrieve the corresponding theta value

       factor_df <- data.frame(features = features, factor = factor)       # Create a data frame with names and retrieved values

       # Extract numeric values using regular expressions
       factor_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_df$features))

       # Sort the data frame based on Iteration and Sigma
       factor_df_sorted <- factor_df[order(factor_df$Index), ]

       # histogram plot
       if (requireNamespace("ggplot2", quietly = TRUE)) {
       plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x = factor, fill = after_stat(x))) +
         geom_histogram(bins=100) +
         labs(x = "Loadings", y = "Frequency", fill = "") + #fill = "Magnitude")
         scale_fill_viridis_c() +  # Use scale_fill_viridis_c() for fill color
         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         theme_minimal()
       } else {
         stop("ggplot2 package is required but not installed.")
       }
       # return
       return(plot)

     } else if (factor_num == 'all'){
       #sim_data <- data.frame(sim_object$concatenated_datasets[[1]])
       #features <- grep("^omic2_feature_", colnames(sim_data), value = TRUE)

       factors <- deltas_names #names(c(sim_object$list_deltas[grep("delta", names(sim_object$list_deltas))]))

       # Initialize plot_list to store plots
       plot_list <- list()

       # Create plots for each factor and store in the plot_list
       for (i in factors) {
         factor_name <- paste("delta", i, sep = "")
         factor <- sim_object$list_deltas[[factor_name]]
         factor_df <- data.frame(features = features, factor = factor)

         # Extract numeric values using regular expressions
         factor_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_df$features))

         # Sort the data frame based on Index
         factor_df_sorted <- factor_df[order(factor_df$Index), ]

         # histogram plot
         if (requireNamespace("ggplot2", quietly = TRUE)) {
         plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x = factor, fill = after_stat(x))) +
           geom_histogram(bins=100) +
           labs(x = "Loadings", y = "Frequency", fill = "") + #fill = "Magnitude")
           scale_fill_viridis_c() +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
           theme_minimal()
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
       message <- paste0("Provide the correct input. Available options: ", paste(deltas_names, collapse = ", "), " and all") # Message
       #print(message) # Print the message and number of available elements
       return(message) # Return the message

     }
   } else if(data == 'integrated'){
     #  # list_betas # list_deltas
     #  # merge the lists
     data <- data.frame(sim_object$concatenated_datasets[[1]])
     features <- colnames(data)

     # Initialize empty vectors to store beta names and numeric values
     list_betas <- sim_object$list_betas
     list_deltas <- sim_object$list_deltas

     # Combine vectors as rows into a data frame
     betas_df <- data.frame(do.call(cbind, list_betas))
     deltas_df <- data.frame(do.call(cbind, list_deltas))

     num_cols_beta <- ncol(betas_df)
     num_cols_delta <- ncol(deltas_df)

     # Extract the last numeric values from the names of the vectors
     betas_names <- as.numeric(gsub("[^0-9]", "", names(sim_object$list_betas)))
     deltas_names <- as.numeric(gsub("[^0-9]", "", names(sim_object$list_deltas)))

     #new_col_beta_names <- paste0("factor", seq_len(num_cols_beta))
     #new_col_delta_names <- paste0("factor", seq_len(num_cols_delta))
     new_col_beta_names <- paste0("factor", betas_names)
     new_col_delta_names <- paste0("factor", deltas_names)
     # Assign new column names
     colnames(betas_df) <- new_col_beta_names
     colnames(deltas_df) <- new_col_delta_names

     # Omic1

     if (factor_num %in% betas_names && factor_num %in% deltas_names) {#if (factor_num %in% seq_along(sim_object$list_betas)) {
       sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

       # Use grep to find column names starting with "omic1_feature_"
       features.omic1 <- grep("^omic1_feature_", colnames(sim_data), value = TRUE)

       factor_omic1 <- sim_object$list_betas[[paste0('beta', factor_num)]] # Retrieve the corresponding theta value

       factor_omic1_df <- data.frame(features = features.omic1, factor = factor_omic1)       # Create a data frame with names and retrieved values

       # Extract numeric values using regular expressions
       factor_omic1_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_omic1_df$features))

       # Sort the data frame based on Iteration and Sigma
       factor_omic1_df_sorted <- factor_omic1_df[order(factor_omic1_df$Index), ]

       # histogram plot
       if (requireNamespace("ggplot2", quietly = TRUE)) {
       plot.a <- ggplot2::ggplot(data = factor_omic1_df_sorted, aes(x = factor, fill = after_stat(x))) +
         geom_histogram(bins=100) +
         labs(x = "Loadings", y = "Frequency", fill = "") + #fill = "Loadings")
         scale_fill_viridis_c("Magnitude") +
         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         theme_bw()+
         ggtitle(paste("Data: omic.one, Factor:", factor_num))
       } else {
         stop("ggplot2 package is required but not installed.")
       }
       sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

       # Use grep to find column names starting with "omic1_feature_"
       features_omic2 <- grep("^omic2_feature_", colnames(sim_data), value = TRUE)

       factor_omic2 <- sim_object$list_deltas[[paste0('delta', factor_num)]] # Retrieve the corresponding theta value

       factor_omic2_df <- data.frame(features = features_omic2, factor = factor_omic2)       # Create a data frame with names and retrieved values

       # Extract numeric values using regular expressions
       factor_omic2_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_omic2_df$features))

       # Sort the data frame based on Iteration and Sigma
       factor_omic2_df_sorted <- factor_omic2_df[order(factor_omic2_df$Index), ]

       # histogram plot
       if (requireNamespace("ggplot2", quietly = TRUE)) {
       plot.b <- ggplot2::ggplot(data = factor_omic2_df_sorted, aes(x = factor, fill = after_stat(x))) +
         geom_histogram(bins=100) +
         labs(x = "Loadings", y = "Frequency", fill = "") + #fill = "Loadings")
         scale_fill_viridis_c("Magnitude") +
         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         theme_bw()+
         ggtitle(paste("Data: omic.two, Factor:", factor_num))
       } else {
         stop("ggplot2 package is required but not installed.")
       }
       # return
       plot_list <- list(plot.a, plot.b)
       if (requireNamespace("gridExtra", quietly = TRUE)) {
       plot_arrange <- gridExtra::grid.arrange(grobs = plot_list, ncol = length(factors))
       } else {
         stop("The 'gridExtra' package is required but not installed. Please install it to use this functionality.")
       }
   }else if(factor_num %in% betas_names && !factor_num %in% deltas_names){
     # Omic1
     sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

     # Use grep to find column names starting with "omic1_feature_"
     features <- grep("^omic1_feature_", colnames(sim_data), value = TRUE)

     if (factor_num %in% betas_names) {#if (factor_num %in% seq_along(sim_object$list_betas)) {

       factor <- sim_object$list_betas[[paste0('beta', factor_num)]] # Retrieve the corresponding theta value

       factor_df <- data.frame(features = features, factor = factor)       # Create a data frame with names and retrieved values

       # Extract numeric values using regular expressions
       factor_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_df$features))

       # Sort the data frame based on Iteration and Sigma
       factor_df_sorted <- factor_df[order(factor_df$Index), ]

       # histogram plot
       if (requireNamespace("ggplot2", quietly = TRUE)) {
       plot <- ggplot2::ggplot(data = factor_df_sorted, aes(x = factor, fill = after_stat(x))) +
         geom_histogram(bins=100) +  # Use geom_bar() instead of geom_histogram()
         labs(x = "Loadings", y = "Frequency", fill = "") + #fill = "Loadings")
         scale_fill_viridis_c() +  # Use scale_fill_viridis_c() for fill color
         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
         theme_minimal()
       } else {
         stop("ggplot2 package is required but not installed.")
       }
       # return
       return(plot)
    }
   }else if(!factor_num %in% betas_names && factor_num %in% deltas_names){
     # Omic2
     # Omic2
     sim_data <- data.frame(sim_object$concatenated_datasets[[1]])

     # Use grep to find column names starting with "omic1_feature_"
     features_omic2 <- grep("^omic2_feature_", colnames(sim_data), value = TRUE)

     factor_omic2 <- sim_object$list_deltas[[paste0('delta', factor_num)]] # Retrieve the corresponding theta value

     factor_omic2_df <- data.frame(features = features_omic2, factor = factor_omic2)       # Create a data frame with names and retrieved values

     # Extract numeric values using regular expressions
     factor_omic2_df$Index <- as.numeric(gsub(".*_(\\d+)$", "\\1", factor_omic2_df$features))

     # Sort the data frame based on Iteration and Sigma
     factor_omic2_df_sorted <- factor_omic2_df[order(factor_omic2_df$Index), ]

     # scatter plot
     if (requireNamespace("ggplot2", quietly = TRUE)) {
     plot.b <- ggplot2::ggplot(data = factor_omic2_df_sorted, aes(x=Index, y=factor, color = factor)) +
       geom_point()+
       labs(x = "Features", y = "Weights", color = paste("")) + #fill = "Loadings")
       theme(
         plot.title = element_text(size = 8, hjust = 0.5)  # Adjust title size and center it
       ) +
       scale_color_viridis_c("Magnitude") +
       theme_bw()+
       ggtitle(paste("Data: Omic.two, Factor:", factor_num))
     } else {
       stop("ggplot2 package is required but not installed.")
     }
     # return
     return(plot.b)

   }
   }
  }
}
# Suppressing the global variable warning
utils::globalVariables(c("Index"))
utils::globalVariables(c("x"))
#################################################################################
# Examples
# plot_weights(output_sim, factor_num = 1, data = 'omic.one', type = 'scatter')
# plot_weights(output_sim, factor_num = 1, data = 'omic.one', type = 'histogram')
# plot_weights(output_sim, factor_num = 3, data = 'omic.one', type = 'scatter')
# plot_weights(output_sim, factor_num = 3, data = 'omic.one', type = 'histogram')
# plot_weights(output_sim, factor_num = 'all', data = 'omic.one', type = 'scatter')
# plot_weights(output_sim, factor_num = 'all', data = 'omic.one', type = 'histogram')

# plot_weights(output_sim, factor_num = 1, data = 'omic.two', type = 'scatter')
# plot_weights(output_sim, factor_num = 1, data = 'omic.two', type = 'histogram')
# plot_weights(output_sim, factor_num = 2, data = 'omic.two', type = 'scatter')
# plot_weights(output_sim, factor_num = 2, data = 'omic.two', type = 'histogram')
# plot_weights(output_sim, factor_num = 'all', data = 'omic.two', type = 'scatter')
# plot_weights(output_sim, factor_num = 'all', data = 'omic.two', type = 'histogram')

# plot_weights(output_sim, factor_num = 3, data = 'integrated', type = 'schttp://127.0.0.1:43149/graphics/31f676fd-25bd-4fdc-929d-1a68f264bf60.pngatter')
# plot_Loadings_hist(output_sim, factor_num = 2, data = 'integrated', type = 'scatter')

