#################################################################
# Function: Configurate the R packages for Parallel-META
# Call: Rscript RM_Config.R
# Authors: Xiaoquan Su
# Last update: 2017-11-01, Xiaoquan Su
# Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
#################################################################

## install necessary libraries
p <- c("Rcpp", "RcppArmadillo", "RcppEigen")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}

invisible(lapply(p, usePackage))

cat("**R Packages Configuration Complete**\n")
