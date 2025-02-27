#' @title Genome-Wide Association Study with GLM and PCA Correction
#' @description Performs a genome-wide association study (GWAS) using a General Linear Model (GLM)
#' with optional Principal Component Analysis (PCA) correction for population structure. Generates
#' p-values, identifies significant SNPs, and produces Manhattan, QQ, and PCA plots if requested.
#' @param pheno A numeric vector or data frame containing phenotype values (n x 1).
#' @param geno A numeric data frame or matrix containing genotype data (n x m), with markers as columns.
#' @param Cov A numeric data frame or matrix containing covariate data (n x t).
#' @param GM A data frame containing SNP information (SNP ID, Chromosome, Position) required for Manhattan plots.
#' @param PCA.M Integer specifying the number of principal components to include as covariates (default: 3).
#' @param cutoff Numeric value for the significance threshold (default: NULL, uses Bonferroni 0.05/m).
#' @param plots Logical indicating whether to generate PCA, Manhattan, and QQ plots (default: FALSE).
#' @param messages Logical indicating whether to print diagnostic messages (default: TRUE).
#' @return A list containing:
#'   \item{P_values}{Numeric vector of p-values for each SNP.}
#'   \item{significant_SNPs}{Integer vector of indices for significant SNPs based on the cutoff.}
#'   \item{PCA}{Matrix of PCA values used in the model, if applicable.}
#' @details This function validates input data, performs PCA to correct for population structure,
#' applies a GLM to calculate p-values, and optionally generates visualizations. It checks for
#' missing values, ensures dimensional consistency, and handles linear dependencies in PCA components
#' using QR decomposition. Plots are skipped if fewer than 3 independent PCs are available.
#' @examples
#' \dontrun{
#' pheno <- rnorm(100)  # Example phenotype
#' geno <- matrix(rbinom(100*1000, 2, 0.3), 100, 1000)  # Example genotype
#' covar <- matrix(rnorm(100*2), 100, 2)  # Example covariates
#' snp_info <- data.frame(SNP = paste0("SNP", 1:1000),
#' Chromosome = sample(1:5, 1000, replace = TRUE),
#' Position = sample(1:1e6, 1000, replace = TRUE))
#' result <- statgenGWAS(pheno = pheno, geno = geno, Cov = covar, GM = snp_info, plots = TRUE)
#' }
#' @export
statgenGWAS <- function(pheno = NULL, geno = NULL, Cov = NULL, GM = NULL,
                        PCA.M = 3, cutoff = NULL, plots = FALSE, messages = TRUE) {
  # Load required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  library(ggplot2)
  library(gridExtra)
  library(data.table)

  # --- Input Validation and Standardization ---
  if (is.null(pheno)) stop("pheno is required.")
  if (is.data.frame(pheno) || is.matrix(pheno)) {
    if (ncol(pheno) != 1) stop("pheno must have 1 column.")
    pheno_vec <- as.numeric(pheno[, 1])
  } else if (is.numeric(pheno)) {
    pheno_vec <- pheno
  } else {
    stop("pheno must be a numeric vector, data frame, or matrix with 1 column.")
  }
  n <- length(pheno_vec)

  if (is.null(geno)) stop("geno is required.")
  if (!is.data.frame(geno) && !is.matrix(geno)) stop("geno must be a data frame or matrix.")
  geno <- as.matrix(geno)
  if (!is.numeric(geno)) stop("geno must be numeric.")
  if (nrow(geno) != n) stop("geno rows must match pheno length.")

  # Store original column names and filter constant/NA-only columns
  original_markers <- colnames(geno)  # Preserve names before filtering
  if (is.null(original_markers)) original_markers <- paste0("SNP", 1:ncol(geno))  # Fallback if no names
  variances <- apply(geno, 2, var, na.rm = TRUE)
  keep_cols <- variances > 0 & !is.na(variances)
  if (sum(keep_cols) < ncol(geno)) {
    if (messages) cat("Removing", ncol(geno) - sum(keep_cols), "constant or NA-only SNPs\n")
    filtered_markers <- original_markers[keep_cols]  # Filtered marker names
    geno <- geno[, keep_cols, drop = FALSE]
    if (!is.null(GM)) GM <- GM[keep_cols, , drop = FALSE]
  } else {
    filtered_markers <- original_markers  # No filtering, keep all names
  }
  m <- ncol(geno)

  if (!is.null(Cov)) {
    if (!is.data.frame(Cov) && !is.matrix(Cov)) stop("Cov must be a data frame or matrix.")
    Cov <- as.matrix(Cov)
    if (!is.numeric(Cov)) stop("Cov must be numeric.")
    if (nrow(Cov) != n) stop("Cov rows must match pheno length.")
    t <- ncol(Cov)
  } else {
    Cov <- matrix(nrow = n, ncol = 0)
    t <- 0
  }

  if (!is.null(GM)) {
    if (!is.data.frame(GM)) stop("GM must be a data frame.")
    if (nrow(GM) != m) stop("GM rows must match geno columns after cleaning.")
    if (!all(c("Chromosome", "Position") %in% colnames(GM))) {
      warning("GM lacks 'Chromosome' or 'Position'; Manhattan plot will be skipped.")
      plots <- FALSE
    }
  } else if (plots) {
    warning("GM is NULL; Manhattan plot requires GM. Skipping plots.")
    plots <- FALSE
  }

  if (any(is.na(pheno_vec)) || any(is.na(geno)) || any(is.na(Cov))) {
    stop("Input data contains missing values. Handle them before running GWAS.")
  }

  # --- PCA Calculation ---
  if (PCA.M > 0) {
    PCA <- prcomp(geno, scale. = TRUE)
    PCA_values_full <- PCA$x[, 1:min(PCA.M, ncol(PCA$x)), drop = FALSE]

    combined <- cbind(Cov, PCA_values_full)
    qr_result <- qr(combined, tol = 1e-10)
    rank <- qr_result$rank
    if (messages) {
      cat("Dimensions of combined:", dim(combined), "\n")
      cat("QR rank:", rank, "\n")
      cat("Pivot indices:", qr_result$pivot, "\n")
    }

    independent_cols <- qr_result$pivot[1:rank]
    pca_cols <- independent_cols[independent_cols > t] - t
    pca_cols <- pca_cols[pca_cols <= PCA.M & pca_cols > 0]

    if (length(pca_cols) > 0) {
      PCA_values <- PCA_values_full[, pca_cols, drop = FALSE]
      if (messages && length(pca_cols) < PCA.M) {
        message(sprintf("Reduced PCA components from %d to %d due to dependence with covariates.",
                        PCA.M, length(pca_cols)))
      }
    } else {
      PCA_values <- NULL
      if (messages) message("No PCA components independent of covariates.")
    }
  } else {
    PCA_values <- NULL
    if (messages) message("PCA.M = 0; skipping PCA.")
  }

  # --- PCA Plots ---
  if (plots && !is.null(PCA_values) && ncol(PCA_values) >= 3) {
    pca_plot_data <- data.frame(PC1 = PCA_values[, 1], PC2 = PCA_values[, 2], PC3 = PCA_values[, 3])
    p1 <- ggplot(pca_plot_data, aes(x = PC1, y = PC2)) + geom_point() + ggtitle("PC1 vs PC2") + theme_minimal()
    p2 <- ggplot(pca_plot_data, aes(x = PC1, y = PC3)) + geom_point() + ggtitle("PC1 vs PC3") + theme_minimal()
    p3 <- ggplot(pca_plot_data, aes(x = PC2, y = PC3)) + geom_point() + ggtitle("PC2 vs PC3") + theme_minimal()
    grid.arrange(p1, p2, p3, nrow = 2, ncol = 2, top = "Principal Component Analysis")
  }

  # --- GLM GWAS: Calculate P-values ---
  P_values <- apply(geno, 2, function(marker) {
    if (max(marker) == min(marker)) return(NA)
    model <- if (is.null(PCA_values)) lm(pheno_vec ~ marker + Cov)
    else lm(pheno_vec ~ marker + Cov + PCA_values)
    summary(model)$coefficients[2, 4]
  })

  # --- Multiple Testing Cutoff ---
  cutoff_final <- if (is.null(cutoff)) 0.05 / m else cutoff
  significant_SNPs <- which(P_values <= cutoff_final)

  # --- Plots (QQ and Manhattan) ---
  if (plots) {
    observed <- -log10(sort(P_values[!is.na(P_values)], decreasing = FALSE))
    expected <- -log10(ppoints(length(observed)))
    qq_plot <- ggplot(data.frame(Expected = expected, Observed = observed),
                      aes(x = Expected, y = Observed)) +
      geom_point() + geom_abline(slope = 1, intercept = 0, color = "red") +
      ggtitle("QQ Plot of GWAS P-values") + theme_minimal()
    print(qq_plot)

    if (!is.null(GM) && all(c("Chromosome", "Position") %in% colnames(GM))) {
      manhattan_data <- data.frame(SNP = filtered_markers, Chromosome = GM$Chromosome,
                                   Position = GM$Position, P_value = P_values)
      manhattan_data <- manhattan_data[order(manhattan_data$Chromosome, manhattan_data$Position), ]
      manhattan_data$logP <- -log10(manhattan_data$P_value)
      manhattan_plot <- ggplot(manhattan_data, aes(x = 1:nrow(manhattan_data), y = logP,
                                                   color = as.factor(Chromosome))) +
        geom_point() + geom_hline(yintercept = -log10(cutoff_final), color = "red", linetype = "dashed") +
        ggtitle("Manhattan Plot") + xlab("SNP Position") + ylab("-log10(P-value)") +
        theme_minimal() + theme(legend.position = "none")
      print(manhattan_plot)
    }
  }

  # --- Output ---
  return(list(P_values = P_values, significant_SNPs = significant_SNPs, PCA = PCA_values,
              Markers = filtered_markers))
}
