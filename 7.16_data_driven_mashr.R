library(mashr)
library(dplyr)

file1 <- read.csv("/Users/cyilin/Desktop/HDRF Data/HDRFData(1)/Akil-1_HPC_GRov - WT.csv")
file2 <- read.csv("/Users/cyilin/Desktop/HDRF Data/HDRFData(1)/Akil-4_DG_EL.GRov - Control.csv")
file3 <- read.csv("/Users/cyilin/Desktop/HDRF Data/HDRFData(1)/Akil-4_DG_LT.GRov - Control.csv")
file4 <- read.csv("/Users/cyilin/Desktop/HDRF Data/HDRFData(1)/Akil-5_HPC_bLR - bHR.csv")

#### 1. Uses hedgesG as effect sizes and varG as standard errors.####
merged_data <- full_join(file1, file2, by = "symbol") %>%
  full_join(., file3, by = "symbol") %>%
  full_join(., file4, by = "symbol")

merged_data <- merged_data %>%
  rowwise() %>%
  mutate(count_non_na = sum(!is.na(c_across(cols = where(is.numeric))))) %>%
  filter(count_non_na >= 3)

merged_data <- na.omit(merged_data)

process_mashr <- function(data) {
  data <- na.omit(data %>% select(starts_with("hedgesG"), starts_with("varG"), symbol))
  betas <- as.matrix(data %>% select(starts_with("hedgesG")))
  std_errors <- as.matrix(data %>% select(starts_with("varG")))
  
  mash_data <- mash_set_data(betas, std_errors)
  U_pca <- cov_pca(mash_data, 2)
  U_ed <- cov_ed(mash_data, U_pca)
  U <- c(U_ed, U_pca)
  m <- mash(mash_data, U, output_lfdr = TRUE)
  
  posterior_means <- m$result$PosteriorMean
  posterior_sds <- m$result$PosteriorSD
  posterior_lfdr <- m$result$lfdr
  posterior_lfsr <- m$result$lfsr
  # Compute p-values from Z-scores
  z_scores <- posterior_means / posterior_sds
  posterior_pvalues <- 2 * pnorm(-abs(z_scores))
  
  
  mashr_results <- data.frame(
    Symbol = data$symbol,
    Posterior_Mean = rowMeans(posterior_means, na.rm = TRUE),
    Posterior_SD = rowMeans(posterior_sds, na.rm = TRUE),
    Posterior_lfdr = rowMeans(posterior_lfdr, na.rm = TRUE),
    Posterior_lfsr = rowMeans(posterior_lfsr, na.rm = TRUE),
    Posterior_PValue = rowMeans(posterior_pvalues, na.rm = TRUE),
    Significant = rowMeans(posterior_pvalues < 0.05, na.rm = TRUE) >= 0.5,
    Significant_Count = rowSums(posterior_pvalues <= 0.05, na.rm = TRUE)
  )

}


mashr_results <- process_mashr(merged_data)
head(mashr_results)

# Keep genes that are significant at least twice among the four estimates
significant_genes_3 <- mashr_results %>% filter(Significant_Count >= 2)
head(significant_genes_3)

#write.csv(significant_genes_3, "/Users/cyilin/Desktop/Sig_genes_hedgeG_data-driven.csv", row.names = FALSE)

# Identify significant genes in meta-analysis
metaresult_1 <- read.csv("/Users/cyilin/Desktop/7.1_revised_4_datasets_mashr_meta_comparison/metaAnalysisResult(1).csv")
meta_significant <- metaresult_1 %>%
  filter(Pvalue <= 0.05) 

comparison_results <- significant_genes_3 %>%
  inner_join(meta_significant, by = "Symbol", suffix = c("_mashr", "_meta"))

pearson_estimate <- cor(comparison_results$Posterior_Mean, comparison_results$Estimate, method = "pearson")
spearman_estimate <- cor(comparison_results$Posterior_Mean, comparison_results$Estimate, method = "spearman")
print(paste("Pearson correlation for estimates:", pearson_estimate))
print(paste("Spearman correlation for estimates:", spearman_estimate))
#"Pearson correlation for estimates: 0.99218755831853"
#"Spearman correlation for estimates: 0.971438902026715"

#### 2. Uses logFC as effect sizes and se as standard errors.####
process_mashr_logFC <- function(data) {
  data <- na.omit(data %>% select(starts_with("logFC"), starts_with("se"), symbol))
  betas <- as.matrix(data %>% select(starts_with("logFC")))
  std_errors <- as.matrix(data %>% select(starts_with("se")))
  
  mash_data <- mash_set_data(betas, std_errors)
  U_pca <- cov_pca(mash_data, 2)
  U_ed <- cov_ed(mash_data, U_pca)
  U <- c(U_ed, U_pca)
  m <- mash(mash_data, U, output_lfdr = TRUE)
  
  posterior_means <- m$result$PosteriorMean
  posterior_sds <- m$result$PosteriorSD
  posterior_lfdr <- m$result$lfdr
  posterior_lfsr <- m$result$lfsr
  # Compute p-values from Z-scores
  z_scores <- posterior_means / posterior_sds
  posterior_pvalues <- 2 * pnorm(-abs(z_scores))
  
  
  mashr_results <- data.frame(
    Symbol = data$symbol,
    Posterior_Mean = rowMeans(posterior_means, na.rm = TRUE),
    Posterior_SD = rowMeans(posterior_sds, na.rm = TRUE),
    Posterior_lfdr = rowMeans(posterior_lfdr, na.rm = TRUE),
    Posterior_lfsr = rowMeans(posterior_lfsr, na.rm = TRUE),
    Posterior_PValue = rowMeans(posterior_pvalues, na.rm = TRUE),
    Significant = rowMeans(posterior_pvalues < 0.05, na.rm = TRUE) >= 0.5,
    Significant_Count = rowSums(posterior_pvalues <= 0.05, na.rm = TRUE)
  )
}

mashr_results_logFC <- process_mashr_logFC(merged_data)
head(mashr_results_logFC)

# Keep genes that are significant at least twice among the four estimates
significant_genes_4 <- mashr_results_logFC %>% filter(Significant_Count >= 2)
head(significant_genes_4)

#write.csv(mashr_results_logFC, "/Users/cathy/Desktop/Sig_genes_logFC_data-driven.csv", row.names = FALSE)

comparison_results_logFC <- significant_genes_4 %>%
  inner_join(meta_significant, by = "Symbol", suffix = c("_mashr", "_meta"))

pearson_estimate_logFC <- cor(comparison_results_logFC$Posterior_Mean, comparison_results_logFC$Estimate, method = "pearson")
spearman_estimate_logFC <- cor(comparison_results_logFC$Posterior_Mean, comparison_results_logFC$Estimate, method = "spearman")
print(paste("Pearson correlation for estimates:", pearson_estimate_logFC))
print(paste("Spearman correlation for estimates:", spearman_estimate_logFC))
#"Pearson correlation for estimates: 0.948720073975938"
#"Spearman correlation for estimates: 0.936990085398678"
