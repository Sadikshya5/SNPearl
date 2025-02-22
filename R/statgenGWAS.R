statgenGWAS <- function(pheno = NULL, geno = NULL, Cov = NULL, GM = NULL, PCA.M = 3, cutoff = NULL, plots = FALSE, messages = TRUE) {
  # Load required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
  if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

  library(ggplot2)
  library(gridExtra)
  library(data.table)

  # --- Question 1: Input Validation for GLM GWAS ---
  if (is.null(pheno) || is.null(geno) || is.null(Cov)) {
    stop("Missing required inputs: pheno, geno, or Cov.")
  }
  if (is.null(GM)) {
    stop("GM (SNP information) is required for Manhattan plot.")
  }

  # Accept pheno as numeric vector or data frame
  if (is.data.frame(pheno)) {
    if (!is.numeric(as.matrix(pheno))) stop("pheno must be numeric.")
    if (ncol(pheno) != 1) stop("pheno must have 1 column.")
    pheno_vec <- pheno[, 1]
  } else if (is.numeric(pheno)) {
    pheno_vec <- pheno
  } else {
    stop("pheno must be a numeric vector or data frame with 1 column.")
  }

  # Validate geno and Cov as data frames or matrices
  if (!is.data.frame(geno) && !is.matrix(geno)) stop("geno must be a data frame or matrix.")
  if (!is.numeric(as.matrix(geno))) stop("geno must be numeric.")
  if (!is.data.frame(Cov) && !is.matrix(Cov)) stop("Cov must be a data frame or matrix.")
  if (!is.numeric(as.matrix(Cov))) stop("Cov must be numeric.")

  # Check dimensions
  n <- length(pheno_vec)
  if (nrow(geno) != n || nrow(Cov) != n) {
    stop("Dimensions mismatch: pheno (n), geno (n x m), Cov (n x t) must have same n.")
  }
  m <- ncol(geno)
  t <- ncol(Cov)

  if (any(is.na(pheno_vec)) || any(is.na(geno)) || any(is.na(Cov))) {
    stop("Input data contains missing values. Handle missing values before running GWAS.")
  }

  # --- Question 2: PCA Calculation and Dependency Exclusion ---
  PCA <- prcomp(geno, scale. = TRUE)
  PCA_values_full <- PCA$x[, 1:min(PCA.M, ncol(PCA$x)), drop = FALSE]

  # Ensure combined is a numeric matrix
  combined <- cbind(Cov, PCA_values_full)
  if (!is.matrix(combined) || !is.numeric(combined)) {
    stop("Combined matrix (Cov + PCA) is not numeric.")
  }

  # QR decomposition with debugging
  qr_result <- qr(combined, tol = 1e-10)
  rank <- qr_result$rank
  if (messages) {
    cat("Dimensions of combined:", dim(combined), "\n")
    cat("QR rank:", rank, "\n")
    cat("Pivot indices:", qr_result$pivot, "\n")
  }

  # Ensure pivot is numeric and subset correctly
  independent_cols <- qr_result$pivot[1:rank]
  if (!is.numeric(independent_cols)) {
    stop("QR pivot indices are not numeric.")
  }

  # Filter PCA columns (indices > t, adjusted to PCA range)
  pca_cols <- independent_cols[independent_cols > t] - t
  pca_cols <- pca_cols[pca_cols <= PCA.M & pca_cols > 0]

  if (length(pca_cols) > 0) {
    PCA_values <- PCA_values_full[, pca_cols, drop = FALSE]
  } else {
    PCA_values <- NULL
    if (messages) message("No PCA components independent of covariates; proceeding without PCs.")
  }

  if (messages && length(pca_cols) < PCA.M) {
    message(sprintf("Reduced PCA components from %d to %d due to linear dependence with covariates.", PCA.M, length(pca_cols)))
  }

  # --- PCA Plots (PC1 vs PC2, PC1 vs PC3, PC2 vs PC3) ---
  if (plots && !is.null(PCA_values) && ncol(PCA_values) >= 3) {
    pca_plot_data <- data.frame(PC1 = PCA_values[, 1], PC2 = PCA_values[, 2], PC3 = PCA_values[, 3])
    p1 <- ggplot(pca_plot_data, aes(x = PC1, y = PC2)) +
      geom_point() + ggtitle("PC1 vs PC2") + theme_minimal()
    p2 <- ggplot(pca_plot_data, aes(x = PC1, y = PC3)) +
      geom_point() + ggtitle("PC1 vs PC3") + theme_minimal()
    p3 <- ggplot(pca_plot_data, aes(x = PC2, y = PC3)) +
      geom_point() + ggtitle("PC2 vs PC3") + theme_minimal()
    grid.arrange(p1, p2, p3, nrow = 2, ncol = 2, top = "Principal Component Analysis")
  } else if (plots && (is.null(PCA_values) || ncol(PCA_values) < 3)) {
    message("Fewer than 3 independent PCs available; PCA plots skipped.")
  }

  # --- GLM GWAS: Calculate P-values ---
  P_values <- apply(geno, 2, function(marker) {
    if (max(marker) == min(marker)) {
      return(NA)
    } else {
      if (is.null(PCA_values)) {
        model <- lm(pheno_vec ~ marker + Cov)
      } else {
        model <- lm(pheno_vec ~ marker + Cov + PCA_values)
      }
      return(summary(model)$coefficients[2, 4])
    }
  })

  # --- Multiple Testing Cutoff ---
  cutoff_final <- if (is.null(cutoff)) 0.05 / m else cutoff
  significant_SNPs <- which(P_values <= cutoff_final)

  # --- QQ Plot ---
  if (plots) {
    observed <- -log10(sort(P_values[!is.na(P_values)], decreasing = FALSE))
    expected <- -log10(ppoints(length(observed)))
    qq_plot <- ggplot(data.frame(Expected = expected, Observed = observed),
                      aes(x = Expected, y = Observed)) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0, color = "red") +
      ggtitle("QQ Plot of GWAS P-values") +
      theme_minimal()
    print(qq_plot)
  }

  # --- Manhattan Plot ---
  if (plots) {
    manhattan_data <- data.frame(
      SNP = colnames(geno),
      Chromosome = GM$Chromosome,
      Position = GM$Position,
      P_value = P_values
    )
    manhattan_data <- manhattan_data[order(manhattan_data$Chromosome, manhattan_data$Position), ]
    manhattan_data$logP <- -log10(manhattan_data$P_value)

    manhattan_plot <- ggplot(manhattan_data, aes(x = 1:nrow(manhattan_data), y = logP, color = as.factor(Chromosome))) +
      geom_point() +
      geom_hline(yintercept = -log10(cutoff_final), color = "red", linetype = "dashed") +
      ggtitle("Manhattan Plot") +
      xlab("SNP Position") +
      ylab("-log10(P-value)") +
      theme_minimal() +
      theme(legend.position = "none")
    print(manhattan_plot)
  }

  # --- Output ---
  return(list(
    P_values = P_values,
    significant_SNPs = significant_SNPs,
    PCA = PCA_values
  ))
}
