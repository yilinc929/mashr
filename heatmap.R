library(mashr)
library(dplyr)
library(tidyr)
library(pheatmap)

file_paths <- c(
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-1_HPC_GRov - WT.csv",
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-4_DG_EL.GRov - Control.csv",
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-4_DG_LT.GRov - Control.csv",
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-5_HPC_bLR - bHR.csv",
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Akil-6_HPC_bLR - bHR.csv",
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/GSE20388_HPC_FSL - FRL.csv",
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/McEwen-4_vDG_FSL - FRL.csv",
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/McEwen-5_vHPC_CORT.WT - Saline.WT.csv",
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Meaney-2_vDG_Control - Stressed.csv",
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Nestler-1_vHPC_Control.Saline - Resilient.Saline.csv",
  "/Users/cathy/Desktop/HDRF Data/HDRFData(1)/Nestler-3_vHPC_Control - Resilient.csv"
)

data_list <- lapply(file_paths, read.csv)

# Assuming the identifier column is 'symbol' or 'ensembl_id'
identifier <- "symbol" 

# Merge datasets by the identifier column, keeping only hedgesG and varG columns
merged_data <- data_list[[1]] %>%
  select(!!sym(identifier), contains("hedgesG"), contains("varG"))

# Extract effect sizes (hedgesG) and variances (varG)
effect_sizes <- merged_data %>%
  select(starts_with("hedgesG")) %>%
  as.matrix()

variances <- merged_data %>%
  select(starts_with("varG")) %>%
  as.matrix()

data <- mash_set_data(Bhat = effect_sizes, Shat = sqrt(variances))

# Estimate data-driven covariances
U.pca <- cov_pca(data, 5)
U.ed <- cov_ed(data, U.pca)
U.canonical <- cov_canonical(data)
U <- c(U.ed, U.canonical)

# Fit mash model
mash_model <- mash(data, U)

# Extract local false sign rate (lfsr)
#lfsr <- get_lfsr(mash_model)

# Convert lfsr to a data frame for easier handling
#lfsr_df <- as.data.frame(lfsr)

# Generate the heatmap
#pheatmap(lfsr_df, 
         #clustering_distance_rows = "euclidean", 
         #clustering_distance_cols = "euclidean", 
         #clustering_method = "average",
         #color = colorRampPalette(c("blue", "white", "red"))(50),
         #main = "Heatmap with Hierarchical Clustering")





####################################################################
###############Display only those following ones####################
symbols_to_display <- c("Dpy19l1", "Ppm1n", "Ntn1", "Slc13a2", "Gram3", "Slc15a2", "Asgr1", 
                        "Scn3a", "Scn1a", "Slc9a1", "Oas3", "Maob", "Sall1", "Zfp408", "Sfi1", 
                        "Pcid1", "Cnih3", "Alox5", "Scn3a1", "Slc39a12", "Mapk11", "Mgst2", 
                        "Ostalpha", "Aldh6a1", "Cplx2", "Fcrls", "Camk4", "Slc12a6", "Kcnh1", 
                        "Dnm1", "Dpysl2", "Ntng2", "Aldh3a1", "Mical2", "Mical3", "Nckap5", 
                        "St6galnac6", "Slc13a3", "Gm5150", "St8sia6")


filtered_data <- merged_data %>%
  filter(!!sym(identifier) %in% symbols_to_display)

effect_sizes <- filtered_data %>%
  select(starts_with("hedgesG")) %>%
  as.matrix()

variances <- filtered_data %>%
  select(starts_with("varG")) %>%
  as.matrix()

data <- mash_set_data(Bhat = effect_sizes, Shat = sqrt(variances))

U.pca <- cov_pca(data, 5)
U.ed <- cov_ed(data, U.pca)
U.canonical <- cov_canonical(data)
U <- c(U.ed, U.canonical)

mash_model <- mash(data, U)

lfsr <- get_lfsr(mash_model)

lfsr_df <- as.data.frame(lfsr)
filenames <- tools::file_path_sans_ext(basename(file_paths))

rownames(lfsr_df) <- filtered_data[[identifier]]
pheatmap(lfsr_df, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "average",
         color = colorRampPalette(c("blue", "white", "yellow"))(50),
         main = "Heatmap with Hierarchical Clustering",
         labels_col = filenames)
