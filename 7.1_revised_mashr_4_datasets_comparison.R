library(mashr)
library(dplyr)

file1 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-1_HPC_GRov - WT.csv")
file2 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-4_DG_EL.GRov - Control.csv")
file3 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-4_DG_LT.GRov - Control.csv")
file4 <- read.csv("/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-5_HPC_bLR - bHR.csv")

merged_data <- full_join(file1, file2, by = "symbol") %>%
  full_join(., file3, by = "symbol") %>%
  full_join(., file4, by = "symbol")

# only symbols that appear in at least three datasets are retained in the final merged dataset.
merged_data <- merged_data %>%
  rowwise() %>%
  mutate(count_non_na = sum(!is.na(c_across(cols = where(is.numeric))))) %>%
  filter(count_non_na >= 3)

merged_data <- na.omit(merged_data)

#### 1. Uses hedgesG as effect sizes and varG as standard errors.####
process_mashr <- function(data) {
  data <- na.omit(data %>% select(starts_with("hedgesG"), starts_with("varG"), symbol))
  betas <- as.matrix(data %>% select(starts_with("hedgesG")))
  std_errors <- as.matrix(data %>% select(starts_with("varG")))
  
  mash_data <- mash_set_data(betas, std_errors)
  U.c <- cov_canonical(mash_data)
  mash_result <- mash(mash_data, U.c, output_lfdr = TRUE)
  
  posterior_means <- mash_result$result$PosteriorMean
  posterior_sds <- mash_result$result$PosteriorSD
  posterior_lfdr <- mash_result$result$lfdr
  posterior_lfsr <- mash_result$result$lfsr
  # Compute p-values from Z-scores
  z_scores <- posterior_means / posterior_sds
  posterior_pvalues <- 2 * pnorm(-abs(z_scores))
  
  
  mashr_results <- data.frame(
    Symbol = data$symbol,
    Posterior_Mean = posterior_means,
    Posterior_SD = posterior_sds,
    Posterior_lfdr = posterior_lfdr,
    Posterior_lfsr = posterior_lfsr,
    Posterior_PValue = posterior_pvalues
  )
  
  #Choose the most significant p-value among the four for each gene
  most_significant_pvalues <- apply(posterior_pvalues, 1, function(row) min(row, na.rm = TRUE))
  mashr_results$Most_Significant_PValue <- most_significant_pvalues
  
  # Compute FDR based on the most significant p-value
  mashr_results$FDR_Most_Significant_PValue <- p.adjust(mashr_results$Most_Significant_PValue, method = "fdr")
  return(mashr_results)
}


mashr_results <- process_mashr(merged_data)
head(mashr_results)

write.csv(mashr_results, "/Users/cathy/Desktop/7.1_revised_mashr_hedgeG_results.csv", row.names = FALSE)

metaresult_1 <- read_csv("/Users/cathy/Desktop/metaAnalysisResult(1).csv")
comparison_results <- mashr_results %>%
  inner_join(metaresult_1, by = "Symbol", suffix = c("_mashr", "_meta"))

# Calculate correlation coefficients
pearson_estimate <- cor(comparison_results$Most_Significant_PValue, comparison_results$Pvalue, method = "pearson")
spearman_estimate <- cor(comparison_results$Most_Significant_PValue, comparison_results$Pvalue, method = "spearman")
print(paste("Pearson correlation for estimates:", pearson_estimate))
print(paste("Spearman correlation for estimates:", spearman_estimate))
#[1] "Pearson correlation for estimates: 0.480334602859262"
#[1] "Spearman correlation for estimates: 0.5239886477232"


#### 2. Uses logFC as effect sizes and se as standard errors.####
process_mashr_logFC <- function(data) {
  data <- na.omit(data %>% select(starts_with("logFC"), starts_with("se"), symbol))
  betas <- as.matrix(data %>% select(starts_with("logFC")))
  std_errors <- as.matrix(data %>% select(starts_with("se")))
  
  mash_data <- mash_set_data(betas, std_errors)
  U.c <- cov_canonical(mash_data)
  mash_result <- mash(mash_data, U.c, output_lfdr = TRUE)
  
  posterior_means <- mash_result$result$PosteriorMean
  posterior_sds <- mash_result$result$PosteriorSD
  posterior_lfdr <- mash_result$result$lfdr
  posterior_lfsr <- mash_result$result$lfsr
  # Compute p-values from Z-scores
  z_scores <- posterior_means / posterior_sds
  posterior_pvalues <- 2 * pnorm(-abs(z_scores))
  
  
  mashr_results <- data.frame(
    Symbol = data$symbol,
    Posterior_Mean = posterior_means,
    Posterior_SD = posterior_sds,
    Posterior_lfdr = posterior_lfdr,
    Posterior_lfsr = posterior_lfsr,
    Posterior_PValue = posterior_pvalues
  )
  
  #Choose the most significant p-value among the four for each gene
  most_significant_pvalues <- apply(posterior_pvalues, 1, function(row) min(row, na.rm = TRUE))
  mashr_results$Most_Significant_PValue <- most_significant_pvalues
  
  # Compute FDR based on the most significant p-value
  mashr_results$FDR_Most_Significant_PValue <- p.adjust(mashr_results$Most_Significant_PValue, method = "fdr")
  return(mashr_results)
}

mashr_results_logFC <- process_mashr_logFC(merged_data)
head(mashr_results_logFC)

write.csv(mashr_results_logFC, "/Users/cathy/Desktop/7.1_revised_mashr_logFC_results.csv", row.names = FALSE)

comparison_results_logFC <- mashr_results_logFC %>%
  inner_join(metaresult_1, by = "Symbol", suffix = c("_mashr", "_meta"))

pearson_estimate_logFC <- cor(comparison_results_logFC$Most_Significant_PValue, comparison_results_logFC$Pvalue, method = "pearson")
spearman_estimate_logFC <- cor(comparison_results_logFC$Most_Significant_PValue, comparison_results_logFC$Pvalue, method = "spearman")
print(paste("Pearson correlation for estimates:", pearson_estimate_logFC))
print(paste("Spearman correlation for estimates:", spearman_estimate_logFC))
#"Pearson correlation for estimates: 0.622340971394345"
#"Spearman correlation for estimates: 0.607141249099528"