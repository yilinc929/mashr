library(mashr)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggpubr)
library(pheatmap)
library(reshape2)

data1 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-1_HPC_GRov - WT.csv")
data2 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-4_DG_EL.GRov - Control.csv")
data3 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-4_DG_LT.GRov - Control.csv")
data4 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-5_HPC_bLR - bHR.csv")
data5 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-6_HPC_bLR - bHR.csv")
data6 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/GSE20388_HPC_FSL - FRL.csv")
data7 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/McEwen-4_vDG_FSL - FRL.csv")
data8 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/McEwen-5_vHPC_CORT.WT - Saline.WT.csv")
data9 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Meaney-2_vDG_Control - Stressed.csv")
data10 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Nestler-1_vHPC_Control.Saline - Resilient.Saline.csv")
data11 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Nestler-3_vHPC_Control - Resilient.csv")


# Define the function to apply mashr to one dataset and return results as a dataframe
process_mashr <- function(data, dataset_name) {
  # Prepare the data
  effects <- as.matrix(data$hedgesG)
  std_errors <- as.matrix(data$varG)
  
  # Apply mashr analysis
  mash_data <- mash_set_data(effects, std_errors)
  U.c <- cov_canonical(mash_data)
  mash_result <- mash(mash_data, U.c)
  
  # Get significant results
  result_indices <- get_significant_results(mash_result)
  significant_genes <- data[result_indices, ]
  
  # Format the results
  formatted_results <- significant_genes %>%
    select(Symbol = symbol, Estimate = hedgesG, SE = varG, Pvalue = pvalue, FDR = fdr)
  
  # Add the dataset identifier if there are any results
  if (nrow(formatted_results) > 0) {
    formatted_results$Dataset <- dataset_name
  }
  
  return(formatted_results)
}

results1 <- process_mashr(data1, "Akil-1_HPC_GRov - WT.csv")

results2 <- process_mashr(data2, "Akil-4_DG_EL.GRov - Control.csv")

results3 <- process_mashr(data3, "Akil-4_DG_LT.GRov - Control.csv")

results4 <- process_mashr(data4, "Akil-5_HPC_bLR - bHR.csv")

results5 <- process_mashr(data5, "Akil-6_HPC_bLR - bHR.csv")

results6 <- process_mashr(data6, "GSE20388_HPC_FSL - FRL.csv")

results7 <- process_mashr(data7, "McEwen-4_vDG_FSL - FRL.csv")

results8 <- process_mashr(data8, "McEwen-5_vHPC_CORT.WT - Saline.WT.csv")

results9 <- process_mashr(data9, "Meaney-2_vDG_Control - Stressed.csv")

results10 <- process_mashr(data10, "Nestler-1_vHPC_Control.Saline - Resilient.Saline.csv")

results11 <- process_mashr(data11, "Nestler-3_vHPC_Control - Resilient.csv")

all_significant_results <- bind_rows(results1, results2, results3, results4, results5, results6, results7, results8, results9, results10, results11)

output_path <- "/Users/cathy/Desktop/results.csv"
write.csv(all_significant_results, output_path, row.names = FALSE)

### Compare the results ###

meta_analysis_results <- read_excel("/Users/cathy/Desktop/HDRF Data/MaleMetaAnalysisResultFinal(1).xlsx")
comparison_results <- merge(all_significant_results, meta_analysis_results, by = "Symbol", suffixes = c("_mashr", "_meta"))

#  Plotting estimates from mashr and meta-analysis
ggplot(comparison_results, aes(x = Estimate_mashr, y = Estimate_meta)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Comparison of mashr and Meta-Analysis Estimates",
       x = "mashr Estimate",
       y = "Meta-Analysis Estimate") +
  theme_minimal()
# Points lying close to the red line suggest a strong agreement between the two sets of estimates. 
# The further the points are from the red line, the less agreement there is between the two methods.

# Pearson and Spearman correlation coefficients
pearson_corr <- cor(comparison_results$Estimate_mashr, comparison_results$Estimate_meta, method = "pearson")
spearman_corr <- cor(comparison_results$Estimate_mashr, comparison_results$Estimate_meta, method = "spearman")

cat("Pearson correlation:", pearson_corr, "\n")
#Pearson correlation: 0.506056 
cat("Spearman correlation:", spearman_corr, "\n")
#Spearman correlation: 0.5058776 
#indicates a moderate positive linear relationship between the mashr and meta-analysis estimates.


# Bland-Altman plot
comparison_results$mean_estimate <- (comparison_results$Estimate_mashr + comparison_results$Estimate_meta) / 2
comparison_results$diff_estimate <- comparison_results$Estimate_mashr - comparison_results$Estimate_meta

ggplot(comparison_results, aes(x = mean_estimate, y = diff_estimate)) +
  geom_point() +
  geom_hline(yintercept = mean(comparison_results$diff_estimate), color = "red", linetype = "dashed") +
  geom_hline(yintercept = mean(comparison_results$diff_estimate) + 1.96 * sd(comparison_results$diff_estimate), color = "blue", linetype = "dashed") +
  geom_hline(yintercept = mean(comparison_results$diff_estimate) - 1.96 * sd(comparison_results$diff_estimate), color = "blue", linetype = "dashed") +
  labs(title = "Bland-Altman Plot",
       x = "Mean of Estimates",
       y = "Difference between Estimates") +
  theme_minimal()
# Red dashed line: the mean difference is close to zero. 
# It suggests that on average, the two methods provide similar estimates.
# Blue dashed lines: the mean difference ± 1.96 times the standard deviation of the differences.
# Visualize the range within which most differences between the two methods fall