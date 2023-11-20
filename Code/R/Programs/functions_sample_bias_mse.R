
# Description:
# This R script includes functions for generating gamma-sar samples, calculating bias and mean squared error (MSE)
# of estimators, and generating plots to visualize the results. The code is designed for statistical analysis
# of several estimators under different sample sizes and replication scenarios.

# Functions:
# - generate_samples: Generates replicated gamma-sar samples with specified parameters.
# - calculate_bias_mse: Computes bias and MSE for estimators across various sample sizes and replications.
# - generate_plot: Generates bias and MSE plots for different mu values and estimators, organized in a grid layout.

# Usage:
# - Adjust parameters as needed and call the functions 

# Dependencies:
# - Ensure that the required packages (e.g., ggplot2) are installed for proper execution.

# Note:
# This code may produce warnings related to deprecated ggplot2 syntax, which can be safely ignored.

# source("../MainFunctions/gamma_sar_sample.r")
# source("../MainFunctions/entropy_gamma_sar.r")


generate_samples <- function(sample_size, replication, mu, L) {
  samples <- vector("list", replication)
  for (r in 1:replication) {
    samples[[r]] <- gamma_sar_sample(L, mu, sample_size)
  }
  return(samples)
}


calculate_variance <- function(sample_sizes, R, B, mu, L, estimators) {
  true_entropy <- entropy_gamma_sar(L, mu)
  
  
  output <- data.frame(SampleSize = integer(0), Estimator = character(0), Variance = numeric(0))
  
  for (ssize in sample_sizes) {
    samples <- generate_samples(ssize, R, mu, L)
    
    for (estimator_name in names(estimators)) {
      estimator <- estimators[[estimator_name]]
      v.entropy <- numeric(R)
      
      for (r in 1:R) {
        sample <- samples[[r]]
        
        if (grepl("Bootstrap", estimator_name)) {
          v.entropy[r] <- estimator(sample, B = B)
        } else {
          v.entropy[r] <- estimator(sample)
        }
      }
      
      variance <- var(v.entropy)
      
      output <- rbind(output, data.frame(SampleSize = ssize, Estimator = estimator_name, Variance = round(variance, 5)))
    }
  }
  
  return(output)
}

calculate_bias_mse <- function(sample_sizes, R, B, mu, L, estimators) {
  true_entropy <- entropy_gamma_sar(L, mu)
  
  output <- data.frame(n = integer(0), Estimator = character(0), Bias = numeric(0), MSE = numeric(0))
  
  for (ssize in sample_sizes) {
    # Generate samples outside the loop
    samples <- generate_samples(ssize, R, mu, L)
    
    for (estimator_name in names(estimators)) {
      estimator <- estimators[[estimator_name]]
      v.entropy <- numeric(R)
      
      for (r in 1:R) {
        sample <- samples[[r]]
        
        if (grepl("Bootstrap", estimator_name)) {
          v.entropy[r] <- estimator(sample, B = B)
        } else {
          v.entropy[r] <- estimator(sample)
        }
      }
      
      mse <- mean((v.entropy - true_entropy)^2)
      bias <- mean(v.entropy) - true_entropy
      
      output <- rbind(output, data.frame(n = ssize, Estimator = estimator_name, Bias = round(bias, 5), MSE = round(mse, 5)))
    }
  }
  
  return(output)
}


# # Function to generate bias and mse
# calculate_bias_mse1 <- function(sample_sizes, R, B, mu, L, estimators) {
#   true_entropy <- entropy_gamma_sar(L, mu)
#   
#   output <- data.frame(mu = numeric(0), n = integer(0), Estimator = character(0), Bias = numeric(0), MSE = numeric(0))
#   
#   for (ssize in sample_sizes) {
#     samples <- generate_samples(ssize, R, mu, L)
#     
#     for (estimator_name in names(estimators)) {
#       estimator <- estimators[[estimator_name]]
#       v.entropy <- numeric(R)
#       
#       for (r in 1:R) {
#         sample <- samples[[r]]
#         
#         if (grepl("Bootstrap", estimator_name)) {
#           v.entropy[r] <- estimator(sample, B = B)
#         } else {
#           v.entropy[r] <- estimator(sample)
#         }
#       }
#       
#       mse <- mean((v.entropy - true_entropy)^2)
#       bias <- mean(v.entropy) - true_entropy
#       
#       output <- rbind(output, data.frame(mu = mu, n = ssize, Estimator = estimator_name, Bias = round(bias, 5), MSE = round(mse, 5)))
#     }
#   }
#   
#   return(output)
# }
# 
# # Function to generate combined table
# generate_combined_table <- function(sample_sizes, R, B, mu_values, L, estimators) {
#   # Lista para almacenar los resultados
#   combined_results <- list()
#   
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse1(sample_sizes, R, B, mu_val, L, estimators)
#     
#     # Agregar los resultados a la lista
#     combined_results <- c(combined_results, results)
#   }
#   
#   # Convertir la lista a un dataframe
#   combined_results_df <- do.call(rbind, combined_results)
#   
#   return(combined_results_df)
# }



generate_plot <- function(sample_sizes, R, B, mu_values, L, estimators, ncol = 2, nrow = 2) {
  # Lista para almacenar los gráficos
  plot_list <- list()
  
  # Bucle para iterar sobre cada valor de mu
  for (mu_val in mu_values) {
    # Calcular resultados para el valor actual de mu
    results <- calculate_bias_mse(sample_sizes, R, B, mu_val, L, estimators)
    #cat("Resultados para mu =", mu_val, "\n")
    
    # Imprimir la tabla de resultados
    #print(results)
    
    df <- as.data.frame(results)
    
    # Plot Bias y MSE en subgráficos separados
    plot_bias <- ggplot(df, aes(x = n, y = Bias, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2) +
      geom_line(linetype = "solid", linewidth = 0.5) +
      labs(y = "Bias", x = "Sample size") +
      guides(color = guide_legend(title = "Estimator")) +
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
    
    # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
    if (mu_val != mu_values[1]) {
      plot_bias = plot_bias + theme(legend.position = "top")
    }
    
    plot_mse <- ggplot(df, aes(x = n, y = MSE, color = Estimator)) +
      geom_hline(yintercept = 0) +
      geom_point(size = 2) +
      geom_line(linetype = "solid", linewidth = 0.5) +
      labs(y = "MSE", x = "Sample size") +
      guides(color = guide_legend(title = "Estimator")) +
      annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
    
    # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
    if (mu_val != mu_values[1]) {
      plot_mse = plot_mse + theme(legend.position = "top")
    }
    
    # Agregar los gráficos a la lista
    plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
  }
  
  # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
  combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
    plot_layout(guides = "collect")
  
  # No mostrar la figura aquí, devolver el objeto combined_plot
  return(combined_plot)
}

# Function to generate tables
generate_table <- function(sample_sizes, R, B, mu_values, L, estimators) {
  # Lista para almacenar las tablas
  table_list <- list()
  
  # Bucle para iterar sobre cada valor de mu
  for (mu_val in mu_values) {
    # Calcular resultados para el valor actual de mu
    results <- calculate_bias_mse(sample_sizes, R, B, mu_val, L, estimators)
    #cat("Resultados para mu =", mu_val, "\n")
    
    # Imprimir la tabla de resultados
    table_list[[as.character(mu_val)]] <- kable(
      results,
      caption = paste("Table for mu =", mu_val),
      format = "latex", 
      booktabs = TRUE
    ) %>%
      kable_styling(latex_options = c("striped", "hold_position"), font_size = 11)
  }
  
  # Return the list of tables
  return(table_list)
}

generate_table2 <- function(sample_sizes, R, B, mu_values, L, estimators) {
  # Crear una tabla para almacenar los resultados de todos los valores de mu
  combined_results <- NULL
  
  # Bucle para iterar sobre cada valor de mu
  for (mu_val in mu_values) {
    # Calcular resultados para el valor actual de mu
    results <- calculate_bias_mse(sample_sizes, R, B, mu_val, L, estimators)
    
    # Agregar una columna con el valor actual de mu en la primera posición
    results <- cbind(mu = mu_val, results)
    
    # Combina los resultados con la tabla principal
    if (is.null(combined_results)) {
      combined_results <- results
    } else {
      combined_results <- rbind(combined_results, results)
    }
  }
  
  # Return the combined table
  return(combined_results)
}



# # Function to generate tables
# generate_table <- function(sample_sizes, R, B, mu_values, L, estimators) {
#   # Lista para almacenar las tablas
#   table_list <- list()
#   
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, B, mu_val, L, estimators)
#     #cat("Resultados para mu =", mu_val, "\n")
#     
#     # Imprimir la tabla de resultados
#     table_list[[as.character(mu_val)]] <- results
#   }
#   
#   # Return the list of tables
#   return(table_list)
# }
# 
# # generate_plot <- function(sample_sizes, R, mu_values, L, estimators, ncol = 2, nrow = 2) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
# 
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, mu_val, L, estimators)
#     cat("Resultados para mu =", mu_val, "\n")
# 
#     # Imprimir la tabla de resultados
#     print(results)
# 
#     df <- as.data.frame(results)
# 
#     # Plot Bias y MSE en subgráficos separados
#     plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "Bias", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
# 
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_bias = plot_bias + theme(legend.position = "top")
#     }
# 
#     plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "MSE", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
# 
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_mse = plot_mse + theme(legend.position = "top")
#     }
# 
#     # Agregar los gráficos a la lista
#     plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
#   }
# 
#   # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
#   combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
#     plot_layout(guides = "collect")
# 
#   # No mostrar la figura aquí, devolver el objeto combined_plot
#   return(combined_plot)
# }
# 



# generate_plot <- function(sample_sizes, R, mu_values, L, estimators, ncol = 2, nrow = 2) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
#   
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, mu_val, L, estimators)
#     cat("Resultados para mu =", mu_val, "\n")
#     
#     # Imprimir la tabla de resultados
#     print(results)
#     
#     df <- as.data.frame(results)
#     
#     # Plot Bias y MSE en subgráficos separados
#     plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "Bias", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = bquote(mu == .(mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#     
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_bias = plot_bias + theme(legend.position = "top")
#     }
#     
#     plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", linewidth = 0.5) +
#       labs(y = "MSE", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = bquote(mu == .(mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#     
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_mse = plot_mse + theme(legend.position = "top")
#     }
#     
#     # Agregar los gráficos a la lista
#     plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
#   }
#   
#   # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
#   combined_plot <- wrap_plots(plot_list, ncol = ncol, nrow = nrow) +
#     plot_layout(guides = "collect")
#   
#   # No mostrar la figura aquí, devolver el objeto combined_plot
#   return(combined_plot)
# }






# generate_plot <- function(sample_sizes, R, mu_values, L, estimators) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
#   
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, mu_val, L, estimators)
#     cat("Resultados para mu =", mu_val, "\n")
#     
#     # Imprimir la tabla de resultados
#     print(results)
#     
#     df <- as.data.frame(results)
#     
#     # Plot Bias y MSE en subgráficos separados
#     plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", size = 0.5) +
#       labs(y = "Bias", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#     
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_bias = plot_bias + theme(legend.position = "top")
#     }
#     
#     plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", size = 0.5) +
#       labs(y = "MSE", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#     
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_mse = plot_mse + theme(legend.position = "top")
#     }
#     
#     # Agregar los gráficos a la lista
#     plot_list[[as.character(mu_val)]] <- list(plot_bias, plot_mse)
#   }
#   
#   # Devolver la lista de gráficos
#   return(plot_list)
# }

# generate_plot <- function(sample_sizes, R, mu_values, L, estimators) {
#   # Lista para almacenar los gráficos
#   plot_list <- list()
# 
#   # Bucle para iterar sobre cada valor de mu
#   for (mu_val in mu_values) {
#     # Calcular resultados para el valor actual de mu
#     results <- calculate_bias_mse(sample_sizes, R, mu_val, L, estimators)
#     cat("Resultados para mu =", mu_val, "\n")
# 
#     # Imprimir la tabla de resultados
#     print(results)
# 
#     df <- as.data.frame(results)
# 
#     # Plot Bias y MSE en subgráficos separados
#     plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", size = 0.5) +
#       labs(y = "Bias", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
# 
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_bias = plot_bias + theme(legend.position = "top")
#     }
# 
#     plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
#       geom_hline(yintercept = 0) +
#       geom_point(size = 2) +
#       geom_line(linetype = "solid", size = 0.5) +
#       labs(y = "MSE", x = "Sample size") +
#       guides(color = guide_legend(title = "Estimator")) +
#       annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
# 
#     # Eliminar la leyenda para todos los gráficos excepto el primero de cada fila
#     if (mu_val != mu_values[1]) {
#       plot_mse = plot_mse + theme(legend.position = "top")
#     }
# 
#     # Agregar los gráficos a la lista
#     plot_list[[as.character(mu_val)]] <- plot_bias + plot_mse
#   }
# 
#   # Organizar los gráficos en una sola figura y añadir leyenda común en la parte inferior
#   combined_plot <- wrap_plots(plot_list, ncol = 2, nrow = 2) +
#     plot_layout(guides = "collect")
# 
#   # Mostrar la figura
#   print(combined_plot)
# }

# ...

# ...

# generate_plots <- function(results, mu_val) {
#   df <- as.data.frame(results)
#   
#   plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
#     geom_hline(yintercept = 0) +
#     geom_point(size = 2) +
#     geom_line(linetype = "solid", size = 0.5) +
#     labs(y = "Bias", x = "Sample size") +
#     guides(color = guide_legend(title = "Estimator")) +
#     annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#   
#   if (mu_val != mu_values[1]) {
#     plot_bias = plot_bias + theme(legend.position = "top")
#   }
#   
#   plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
#     geom_hline(yintercept = 0) +
#     geom_point(size = 2) +
#     geom_line(linetype = "solid", size = 0.5) +
#     labs(y = "MSE", x = "Sample size") +
#     guides(color = guide_legend(title = "Estimator")) +
#     annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("mu == %s", mu_val)), hjust = 1.08, vjust = 1.3, size = 3)
#   
#   if (mu_val != mu_values[1]) {
#     plot_mse = plot_mse + theme(legend.position = "top")
#   }
#   
#   return(plot_bias + plot_mse)
# }
# 
# 
# calculate_bias_mse <- function(sample_sizes, R, mu, L, estimators) {
#   true_entropy <- entropy_gamma_sar(L, mu)
#   
#   output <- data.frame(SampleSize = integer(0), Estimator = character(0), Bias = numeric(0), MSE = numeric(0))
#   
#   for (ssize in sample_sizes) {
#     # Generate samples outside the loop
#     samples <- generate_samples(ssize, R, mu, L)
#     
#     for (estimator_name in names(estimators)) {
#       estimator <- estimators[[estimator_name]]
#       v.entropy <- numeric(R)
#       
#       for (r in 1:R) {
#         sample <- samples[[r]]
#         
#         if (grepl("Bootstrap", estimator_name)) {
#           v.entropy[r] <- estimator(sample, B = 10)
#         } else {
#           v.entropy[r] <- estimator(sample)
#         }
#       }
#       
#       mse <- mean((v.entropy - true_entropy)^2)
#       bias <- mean(v.entropy) - true_entropy
#       
#       output <- rbind(output, data.frame(SampleSize = ssize, Estimator = estimator_name, Bias = round(bias, 5), MSE = round(mse, 5)))
#     }
#   }
#   
#   return(output)
# }

