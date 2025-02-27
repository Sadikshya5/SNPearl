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
# Create data frame with filtered markers and P-values
result_df <- data.frame(Marker = result$Markers, P_value = result$P_values)
cat("Result_df dimensions:", dim(result_df), "\n")  # Should be 2953 x 2

# Save to Excel
write_xlsx(list("P_Values" = result_df), "statgenGWAS_Results.xlsx")
cat("Results saved to statgenGWAS_Results.xlsx\n")

