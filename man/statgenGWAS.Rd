% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statgenGWAS.R
\name{statgenGWAS}
\alias{statgenGWAS}
\title{Association Analysis with PCA, Manhattan Plot, and QQ Plot}
\usage{
statgenGWAS(
  pheno = NULL,
  geno = NULL,
  Cov = NULL,
  GM = NULL,
  PCA.M = 3,
  cutoff = NULL,
  plots = FALSE,
  messages = TRUE
)
}
\arguments{
\item{pheno}{A numeric vector of phenotype values.}

\item{geno}{A numeric matrix of genotype data, with markers as columns.}

\item{Cov}{A numeric matrix of covariate data.}

\item{GM}{A data frame containing SNP information (SNP ID, chromosome, position).}

\item{PCA.M}{An integer indicating the number of principal components to include as cofactors (default is 3).}

\item{cutoff}{A numeric value for the significance threshold (optional).}

\item{plots}{A logical value indicating whether to generate PCA, Manhattan, and QQ plots (default is FALSE).}

\item{messages}{A logical value indicating whether to print messages during the execution (default is FALSE).}

\item{QTN.position}{A vector indicating the positions of known quantitative trait nucleotides (optional).}

\item{print}{A logical value indicating whether to print the results (default is FALSE).}

\item{trait}{A character string indicating the name of the trait (default is "unknown").}
}
\value{
A list containing p-values, significant SNPs, and PCA results.
}
\description{
This function performs genome-wide association studies (GWAS) with the option to include PCA components,
generates Manhattan and QQ plots, and identifies significant SNPs.
}
