# Load required libraries
library(ggplot2)
library(data.table)
library(pROC)
source("http://zzlab.net/StaGen/2020/R/G2P.R")
source("http://zzlab.net/StaGen/2020/R/GWASbyCor.R")

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

# GWASbyCor function (as provided)
GWASbyCor <- function(X, y) {
  n <- nrow(X)
  r <- cor(y, X)
  t <- r / sqrt((1 - r^2) / (n - 2))
  p <- 2 * (1 - pt(abs(t), n - 2))
  zeros <- p == 0
  p[zeros] <- 1e-10
  return(p)
}

# Simulation function with per-replicate stats for t-test
simulate_GWAS_comparison <- function(geno, Cov, GM, n_qtn = 10, effect_size = 0.5, n_reps = 30) {
  set.seed(123)
  n <- nrow(geno)
  m <- ncol(geno)
  thresholds <- seq(0, 0.1, by = 0.001)
  results <- list(statgenGWAS = list(power = numeric(length(thresholds)), fdr = numeric(length(thresholds)), time = numeric(n_reps),
                                     power_rep = numeric(n_reps), fdr_rep = numeric(n_reps)),
                  GWASbyCor = list(power = numeric(length(thresholds)), fdr = numeric(length(thresholds)), time = numeric(n_reps),
                                   power_rep = numeric(n_reps), fdr_rep = numeric(n_reps)))

  for (rep in 1:n_reps) {
    # Select QTNs from your genotype data
    qtn_cols <- sample(1:m, n_qtn)
    qtn_effect <- rnorm(n_qtn, mean = effect_size, sd = 0.1)

    # Simulate phenotype
    pheno <- rowSums(t(t(geno[, qtn_cols]) * qtn_effect)) + rowSums(Cov) + rnorm(n, sd = 1)

    # Run statgenGWAS
    start_time <- Sys.time()
    stat_res <- statgenGWAS(pheno = pheno, geno = geno, Cov = Cov, GM = GM)
    results$statgenGWAS$time[rep] <- as.numeric(Sys.time() - start_time)

    # Run GWASbyCor
    start_time <- Sys.time()
    cor_pvals <- GWASbyCor(geno, pheno)
    results$GWASbyCor$time[rep] <- as.numeric(Sys.time() - start_time)

    # Calculate power and FDR per replicate (using a fixed threshold, e.g., 0.05)
    true_positives <- qtn_cols
    thresh <- 0.05  # Fixed threshold for per-replicate stats

    # statgenGWAS
    sig_stat <- which(stat_res$P_values <= thresh)
    tp_stat <- length(intersect(sig_stat, true_positives))
    fp_stat <- length(setdiff(sig_stat, true_positives))
    results$statgenGWAS$power_rep[rep] <- tp_stat / n_qtn
    results$statgenGWAS$fdr_rep[rep] <- if (length(sig_stat) > 0) fp_stat / length(sig_stat) else 0

    # GWASbyCor
    sig_cor <- which(cor_pvals <= thresh)
    tp_cor <- length(intersect(sig_cor, true_positives))
    fp_cor <- length(setdiff(sig_cor, true_positives))
    results$GWASbyCor$power_rep[rep] <- tp_cor / n_qtn
    results$GWASbyCor$fdr_rep[rep] <- if (length(sig_cor) > 0) fp_cor / length(sig_cor) else 0

    # ROC curve data (averaged across replicates)
    for (i in seq_along(thresholds)) {
      thresh <- thresholds[i]

      # statgenGWAS
      sig_stat <- which(stat_res$P_values <= thresh)
      tp_stat <- length(intersect(sig_stat, true_positives))
      fp_stat <- length(setdiff(sig_stat, true_positives))
      power_stat <- tp_stat / n_qtn
      fdr_stat <- if (length(sig_stat) > 0) fp_stat / length(sig_stat) else 0
      results$statgenGWAS$power[i] <- results$statgenGWAS$power[i] + power_stat / n_reps
      results$statgenGWAS$fdr[i] <- results$statgenGWAS$fdr[i] + fdr_stat / n_reps

      # GWASbyCor
      sig_cor <- which(cor_pvals <= thresh)
      tp_cor <- length(intersect(sig_cor, true_positives))
      fp_cor <- length(setdiff(sig_cor, true_positives))
      power_cor <- tp_cor / n_qtn
      fdr_cor <- if (length(sig_cor) > 0) fp_cor / length(sig_cor) else 0
      results$GWASbyCor$power[i] <- results$GWASbyCor$power[i] + power_cor / n_reps
      results$GWASbyCor$fdr[i] <- results$GWASbyCor$fdr[i] + fdr_cor / n_reps
    }
  }

  return(results)
}

# Run simulation with your data
sim_results <- simulate_GWAS_comparison(geno = X, Cov = C, GM = GM, effect_size = 0.5, n_reps = 30)

# Plot ROC curve (Power vs FDR)
roc_data <- data.frame(
  FDR = c(sim_results$statgenGWAS$fdr, sim_results$GWASbyCor$fdr),
  Power = c(sim_results$statgenGWAS$power, sim_results$GWASbyCor$power),
  Method = rep(c("statgenGWAS", "GWASbyCor"), each = length(sim_results$statgenGWAS$fdr))
)
ggplot(roc_data, aes(x = FDR, y = Power, color = Method)) +
  geom_line() +
  labs(title = "ROC Curve: Power vs FDR (Using Own Data)", x = "False Discovery Rate", y = "Power (TPR)") +
  theme_minimal()


# Summarize AUC
auc_stat <- sum(diff(sim_results$statgenGWAS$fdr) * sim_results$statgenGWAS$power[-1])
auc_cor <- sum(diff(sim_results$GWASbyCor$fdr) * sim_results$GWASbyCor$power[-1])
cat("AUC (Power vs FDR) - statgenGWAS:", auc_stat, "GWASbyCor:", auc_cor, "\n")
cat("Mean Time - statgenGWAS:", mean(sim_results$statgenGWAS$time), "GWASbyCor:", mean(sim_results$GWASbyCor$time), "\n")

# Perform t-tests
# Power
power_ttest <- t.test(sim_results$statgenGWAS$power_rep, sim_results$GWASbyCor$power_rep, paired = TRUE, alternative = "greater")
cat("Paired t-test for Power (statgenGWAS > GWASbyCor):\n")
print(power_ttest)

# FDR
fdr_ttest <- t.test(sim_results$statgenGWAS$fdr_rep, sim_results$GWASbyCor$fdr_rep, paired = TRUE, alternative = "less")
cat("Paired t-test for FDR (statgenGWAS < GWASbyCor):\n")
print(fdr_ttest)

