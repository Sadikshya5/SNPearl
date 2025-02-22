library(devtools)  # Install devtools if not installed
devtools::install("C:/Users/sadik/OneDrive - Washington State University (email.wsu.edu)/Course/Crop545_Statgenomics/Assignment 1/SNPearl")
library(SNPearl)
# Your data loading code
phenotype <- read.table("phenotype.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
genotype <- read.table("GAPIT.genotype.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
covariates <- read.table("Barley_covariates.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
GM <- read.table("GM.txt", header = TRUE)

cat("Phenotype:", dim(phenotype), "Genotype:", dim(genotype), "Covariates:", dim(covariates), "\n")

y <- as.numeric(as.matrix(phenotype[, -1, drop = FALSE]))
X <- as.matrix(sapply(genotype[, -1], as.numeric))
C <- as.matrix(sapply(covariates[, -1], as.numeric))

str(y); str(X); str(C)
cat("Missing values - y:", sum(is.na(y)), "X:", sum(is.na(X)), "C:", sum(is.na(C)), "\n")

# Run the function with debugging output
result <- statgenGWAS(pheno = y, geno = X, Cov = C, GM = GM, plots = TRUE, messages = TRUE)


#Using statgenwithoutPCA to get p_values
#p_values<- statgenGWAS_withoutpca (y, X, C, PCA.M = 3, cutoff = NULL, plots = FALSE, verbose = TRUE)
#head(p_values)


###################
# Load necessary libraries to download excel file with the p- values
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
library(writexl)
# Convert to a data frame
#p_values_df <- data.frame(Marker = colnames(X), P_value = p_values_result$P_values)
# Write to Excel
#write_xlsx(list("P_Values" = p_values_df), "GWAS_p_values.xlsx")
#
# Extract p-values
result_df <- data.frame(Marker = colnames(X), P_value = result$P_values)
# Save to Excel
write_xlsx(list("P_Values" = result_df), "statgenGWAS_Results.xlsx")





