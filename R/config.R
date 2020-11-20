#' install necessary libraries
#' @param p the list of necessary packages
#' @param usePackage

p <- c("Rcpp", "RcppArmadillo", "RcppEigen")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}

invisible(lapply(p, usePackage))

cat("**R Packages Configuration Complete**\n")
