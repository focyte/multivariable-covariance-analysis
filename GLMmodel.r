# Load necessary libraries
library(ggplot2)
library(dplyr)
library(caret)
library(readr)
library(broom)
library(boot)
library(relaimpo)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##                            ███╗   ███╗ ██████╗ █████╗
##                            ████╗ ████║██╔════╝██╔══██╗
##                            ██╔████╔██║██║     ███████║
##                            ██║╚██╔╝██║██║     ██╔══██║
##                            ██║ ╚═╝ ██║╚██████╗██║  ██║
##                            ╚═╝     ╚═╝ ╚═════╝╚═╝  ╚═╝

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read the data from TSV files in the results folder
x <- read_tsv("x.tsv")
y <- read_tsv("y.tsv")

# Combine the data into one dataframe
data <- cbind(x, y)

# Define the pre-processing steps
preprocess_params <-
  preProcess(data[, -ncol(data)],
  method = c("center", "scale", "YeoJohnson", "nzv"))

# Normalize the data
normalized_data <- predict(preprocess_params, newdata = data[, -ncol(data)])

# Combine the normalized predictors with the response variable
final_data <- cbind(normalized_data, BuffaScore = data$BuffaScore)

# Fit a GLM model
glm_model <- glm(BuffaScore ~ ., data = final_data)

# Save the summary of the GLM model to a text file in the results folder
summary_text <- capture.output(summary(glm_model))
writeLines(summary_text, paste0("glm_model_summary_hypoxia.txt"))

# Extract coefficients and convert to a tidy dataframe
coefficients <- tidy(glm_model)

# Save the coefficients to a CSV file in the results folder
write_csv(coefficients, paste0("glm_coefficients_hypoxia.csv"))

# Print the coefficients to the console
print(coefficients)

# Function to extract coefficients and plot, saving to a PDF
get_coef <- function(glm_model) {
  # Extract coefficients and p-values
  coe <- coef(glm_model)
  pvalue <- coef(summary(glm_model))[, 4]

  # Remove intercept and sort coefficients
  coe <- coe[-1]
  coe <- sort(coe)

  # Create a data frame with coefficients and p-values
  coe_df <- data.frame(g = names(coe), coe = coe, pval = pvalue[names(coe)])
  coe_df$pval[coe_df$pval > 0.05] <- 1

  # Convert 'g' to a factor with sorted levels
  coe_df$g <- factor(coe_df$g, levels = coe_df$g)

  # Plot using ggplot
  plot <- ggplot(coe_df, aes(x = g, y = coe, fill = coe, size = pval < 0.05)) +
    geom_bar(stat = 'identity', color = 'black') +
    theme_bw() +
    coord_flip() +
    scale_fill_gradient2(midpoint = 0, 
      low = "blue", mid = "white", high = "red") +
    scale_size_manual(values = c(0, 1))

  # Save the plot as a PDF in the results folder
  ggsave(filename = paste0("coefficients_plot.pdf"),
    plot = plot, device = "pdf")

  # Print the plot
  print(plot)

  # Return the data frame of coefficients
  return(coe_df)
}

# Call the function with your glm_model
coefficients <- get_coef(glm_model)

# Set seed for reproducibility
set.seed(1910)

# Function to perform bootstrapping analysis for 
# variable importance and calculate %R2 contribution
perform_bootstrapping_relimp <- function(glm_model, b = 1000) {
  # Define bootstrapping function to calculate relative importance
  boot_fun <- function(data, indices) {
    fit <- tryCatch(glm(BuffaScore ~ ., 
      data = data[indices, ]), error = function(e) return(NULL))
    if (!is.null(fit)) {
      relimp <- tryCatch(calc.relimp(fit, type = "lmg")
        $lmg, error = function(e) return(rep(NA, ncol(data) - 1)))
    } else {
      relimp <- rep(NA, ncol(data) - 1)
    }
    return(relimp)
  }
  # Perform bootstrapping
  boot_results <- boot(data = final_data, statistic = boot_fun, R = b)

  # Remove any rows with NA values to prevent issues in summarization
  boot_results_clean <- boot_results$t[complete.cases(boot_results$t), ]

  if (nrow(boot_results_clean) == 0) {
    stop("All bootstrapping attempts failed. Please check the model and data.")
  }

  # Extract mean and standard deviation of relative importance
  relimp_mean <- apply(boot_results_clean, 2, mean, na.rm = TRUE)
  relimp_sd <- apply(boot_results_clean, 2, sd, na.rm = TRUE)

  # Calculate %R2 contribution for each predictor
  # Sum of all contributions gives the total R2
  total_r2 <- sum(relimp_mean)  
  # Calculate percentage contribution
  percent_r2_contribution <- (relimp_mean / total_r2) * 100  

  # Calculate the standard deviation of %R2 contribution for each predictor
  percent_r2_contribution_sd <- (1 / total_r2) * 100 * relimp_sd

  # Combine results into a data frame with the new standard deviation
  relimp_summary <- data.frame(
    # Get predictor names (excluding response variable)
    term = colnames(final_data)[-ncol(final_data)],
    Mean_Contribution = relimp_mean,
    Std_Dev_Contribution = relimp_sd,
    Percent_R2_Contribution = percent_r2_contribution,
    Std_Dev_Percent_R2_Contribution = percent_r2_contribution_sd
  )

  # Save the summary results to a CSV file
  output_csv <- paste0("relimp_summary_hypoxia.csv")
  write.csv(relimp_summary, output_csv, row.names = FALSE)

  # Return the summary data frame
  return(relimp_summary)
}

# Call the function with your glm_model to perform bootstrapping
relimp_results <- perform_bootstrapping_relimp(glm_model, b = 1000)

# Create a combined CSV file from the two result files
# Join the two dataframes by term and select the required columns
coefficients <- coefficients %>%
  rename(
    term = g,
    Beta = coe,
    )

combined_results <- merge(x=relimp_results,y=coefficients,
  by.x="term",by.y = "term",all= T)
print(combined_results)
combined_results <- combined_results %>%
  rename(
    Gene = term,
    Coef_Mean = Mean_Contribution,
    Coef_SD = Std_Dev_Contribution,
    R2 = Percent_R2_Contribution,
    R2_SD = Std_Dev_Percent_R2_Contribution,
  )

print(combined_results)

# Write the combined results to a new CSV file
write_csv(combined_results, "results.csv")

# Modify R2 values based on Beta: if Beta is negative, make R2 negative
combined_results$R2 <- ifelse(combined_results$Beta 
  < 0, -abs(combined_results$R2), abs(combined_results$R2))

# Create the plot
plot <- ggplot(combined_results, aes(x = Gene, y = R2, fill = Beta)) +
  geom_bar(stat = "identity", color = "black") +  # Add bars for the values
  theme_bw() +
  coord_flip() +
  scale_fill_gradient(low = "#0026ffff", high = "red") + 
  labs(title = "Significant Coefficients", 
    fill = "Standardised β Coefficient") +
  scale_y_continuous(breaks=seq(-50, 100, 25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot to a PDF file
ggsave(filename = "R2_plot.pdf", plot = plot, device = "pdf")

# Print the plot
print(plot)
