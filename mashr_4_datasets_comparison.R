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

#### 1. Uses hedgesG as effect sizes and varG as standard errors.####

process_mashr <- function(data) {
  betas <- as.matrix(data %>% select(starts_with("hedgesG")))
  std_errors <- as.matrix(data %>% select(starts_with("varG")))
  
  mash_data <- mash_set_data(betas, std_errors)
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

mashr_results <- process_mashr(merged_data)
head(mashr_results)

# Summarize estimates and SEs after the function
weights <- 1 / (mashr_results[, c("SE1", "SE2", "SE3", "SE4")]^2)
combined_estimate <- rowSums(mashr_results[, c("Estimate1", "Estimate2", "Estimate3", "Estimate4")] * weights) / rowSums(weights)
combined_se <- sqrt(1 / rowSums(weights))

# Calculate p-values from the combined estimates and SEs
combined_pvalue <- 2 * pnorm(-abs(combined_estimate / combined_se))

final_results_1 <- data.frame(
  Symbol = mashr_results$Symbol,
  Estimate = combined_estimate,
  SE = combined_se,
  Pvalue = combined_pvalue,
  FDR = p.adjust(combined_pvalue, method = "fdr")
)

head(final_results_1)

write.csv(final_results_1, "/Users/cathy/Desktop/mashr_results_hedgesG.csv", row.names = FALSE)

metaresult_1 <- read_csv("/Users/cathy/Desktop/metaAnalysisResult.csv")
comparison_results <- final_results_1 %>%
  inner_join(metaresult_1, by = "Symbol", suffix = c("_mashr", "_meta"))

# Calculate correlation coefficients
pearson_estimate <- cor(comparison_results$Estimate_mashr, comparison_results$Estimate_meta, method = "pearson")
spearman_estimate <- cor(comparison_results$Estimate_mashr, comparison_results$Estimate_meta, method = "spearman")
print(paste("Pearson correlation for estimates:", pearson_estimate))
print(paste("Spearman correlation for estimates:", spearman_estimate))
#"Pearson correlation for estimates: 0.815612934549745"
#"Spearman correlation for estimates: 0.804221391656552"

library(ggplot2)
library(ggpubr)
ggplot(comparison_results, aes(x = Estimate_mashr, y = Estimate_meta)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", col = "red") +
  ggtitle("Scatter Plot of Estimates") +
  xlab("mashr Estimates") +
  ylab("Meta-Analysis Estimates")

ggplot(comparison_results, aes(x = SE_mashr, y = SE_meta)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", col = "red") +
  ggtitle("Scatter Plot of Standard Errors") +
  xlab("mashr SE") +
  ylab("Meta-Analysis SE")

#### 2. Uses logFC as effect sizes and se as standard errors.####
process_mashr_logFC <- function(data) {
  betas <- as.matrix(data %>% select(starts_with("logFC")))
  std_errors <- as.matrix(data %>% select(starts_with("se")))
  
  mash_data <- mash_set_data(betas, std_errors)
  U.c <- cov_canonical(mash_data)
  mash_result <- mash(mash_data, U.c)
  
  result_indices <- get_significant_results(mash_result)
  if (length(result_indices) == 0) {
    return(data.frame()) 
  }
  significant_genes <- data[result_indices, ]
  
  formatted_results <- significant_genes %>%
    select(Symbol = symbol,
           Estimate = starts_with("hedgesG"),
           SE = starts_with("varG"))
  return(formatted_results)
}

mashr_results_logFC <- process_mashr_logFC(merged_data)
head(mashr_results_logFC)

weights <- 1 / (mashr_results_logFC[, c("SE1", "SE2", "SE3", "SE4")]^2)
combined_estimate <- rowSums(mashr_results_logFC[, c("Estimate1", "Estimate2", "Estimate3", "Estimate4")] * weights) / rowSums(weights)
combined_se <- sqrt(1 / rowSums(weights))

combined_pvalue <- 2 * pnorm(-abs(combined_estimate / combined_se))

final_results_logFC <- data.frame(
  Symbol = mashr_results_logFC$Symbol,
  Estimate = combined_estimate,
  SE = combined_se,
  Pvalue = combined_pvalue,
  FDR = p.adjust(combined_pvalue, method = "fdr")
)

head(final_results_logFC)

write.csv(final_results_logFC, "/Users/cathy/Desktop/mashr_results_logFC.csv", row.names = FALSE)

comparison_results_logFC <- final_results_logFC %>%
  inner_join(metaresult_1, by = "Symbol", suffix = c("_mashr", "_meta"))

pearson_estimate_logFC <- cor(comparison_results_logFC$Estimate_mashr, comparison_results_logFC$Estimate_meta, method = "pearson")
spearman_estimate_logFC <- cor(comparison_results_logFC$Estimate_mashr, comparison_results_logFC$Estimate_meta, method = "spearman")
print(paste("Pearson correlation for estimates:", pearson_estimate_logFC))
print(paste("Spearman correlation for estimates:", spearman_estimate_logFC))
#"Pearson correlation for estimates: 0.864207076070862"
#"Spearman correlation for estimates: 0.801270326023692"

ggplot(comparison_results_logFC, aes(x = Estimate_mashr, y = Estimate_meta)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", col = "red") +
  ggtitle("Scatter Plot of Estimates (logFC)") +
  xlab("mashr Estimates") +
  ylab("Meta-Analysis Estimates")

ggplot(comparison_results_logFC, aes(x = SE_mashr, y = SE_meta)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", col = "red") +
  ggtitle("Scatter Plot of Standard Errors (logFC)") +
  xlab("mashr SE") +
  ylab("Meta-Analysis SE")