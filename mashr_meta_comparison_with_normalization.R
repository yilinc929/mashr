library(mashr)
library(dplyr)

file1 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-1_HPC_GRov - WT.csv")
file2 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-4_DG_EL.GRov - Control.csv")
file3 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-4_DG_LT.GRov - Control.csv")
file4 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-5_HPC_bLR - bHR.csv")

merged_data <- file1 %>%
  inner_join(file2, by = "symbol") %>%
  inner_join(file3, by = "symbol") %>%
  inner_join(file4, by = "symbol")

normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

remove_outliers <- function(data, columns, threshold = 1.5) {
  for (column in columns) {
    Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
    Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - threshold * IQR
    upper_bound <- Q3 + threshold * IQR
    data <- data %>% filter(data[[column]] >= lower_bound & data[[column]] <= upper_bound)
  }
  return(data)
}

columns_to_process <- c("hedgesG.x", "hedgesG.y", "varG.x", "varG.y")

processed_data <- merged_data %>%
  mutate(across(all_of(columns_to_process), normalize)) %>%
  remove_outliers(columns_to_process)

processed_data <- processed_data %>%
  mutate(across(starts_with("varG"), ~ ifelse(. < .Machine$double.eps, .Machine$double.eps, .)))

head(processed_data)

#### 1. Uses hedgesG as effect sizes and varG as standard errors.####

process_mashr_1 <- function(data) {
  betas <- as.matrix(data %>% select(starts_with("hedgesG")))
  std_errors <- as.matrix(data %>% select(starts_with("varG")))
  
  mash_data <- mash_set_data(betas, std_errors, zero_Shat_reset = .Machine$double.eps)
  U.c <- cov_canonical(mash_data)
  mash_result <- mash(mash_data, U.c)
  
  # Get significant results
  result_indices <- get_significant_results(mash_result)
  if (length(result_indices) == 0) {
    return(data.frame())  # Return an empty dataframe if no significant results
  }
  significant_genes <- data[result_indices, ]
  
  formatted_results <- significant_genes %>%
    select(Symbol = symbol,
           Estimate = starts_with("hedgesG"),
           SE = starts_with("varG"))
  return(formatted_results)
}

mashr_results <- process_mashr_1(processed_data)
head(mashr_results)

# Summarize estimates and SEs after the function
weights <- 1 / (mashr_results[, c("SE1", "SE2", "SE3", "SE4")]^2)
combined_estimate <- rowSums(mashr_results[, c("Estimate1", "Estimate2", "Estimate3", "Estimate4")] * weights) / rowSums(weights)
combined_se <- sqrt(1 / rowSums(weights))

# Calculate p-values from the combined estimates and SEs
combined_pvalue <- 2 * pnorm(-abs(combined_estimate / combined_se))

final_results_1_1 <- data.frame(
  Symbol = mashr_results$Symbol,
  Estimate = combined_estimate,
  SE = combined_se,
  Pvalue = combined_pvalue,
  FDR = p.adjust(combined_pvalue, method = "fdr")
)

head(final_results_1_1)

#write.csv(final_results_1_1, "/Users/cathy/Desktop/mashr_results_hedgesG.csv", row.names = FALSE)

#metaresult_1 <- read_csv("/Users/cathy/Desktop/metaAnalysisResult.csv")
comparison_results <- final_results_1_1 %>%
  inner_join(metaresult_1, by = "Symbol", suffix = c("_mashr", "_meta"))

# Calculate correlation coefficients
pearson_estimate <- cor(comparison_results$Estimate_mashr, comparison_results$Estimate_meta, method = "pearson")
spearman_estimate <- cor(comparison_results$Estimate_mashr, comparison_results$Estimate_meta, method = "spearman")
print(paste("Pearson correlation for estimates:", pearson_estimate))
print(paste("Spearman correlation for estimates:", spearman_estimate))
#"Pearson correlation for estimates: 0.815612934549745"
#"Spearman correlation for estimates: 0.804221391656552"