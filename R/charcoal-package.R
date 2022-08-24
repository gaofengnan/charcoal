#' charcoal: An R package implementing the charcoal algorithms for localizing
#' changepoints in high-dimensional linear regression (via complementary
#' sketching)
#' 
#'
#' @docType package
#' @name charcoal
#'
#' @description The charcoal package collects the implementation for the novel 
#' methodology #' 'charcoal' for changepoint localization in high-dimensional 
#' linear regression, introduced in Gao and Wang (2022).
#' @references Gao, F. and Wang, T. (2022) Sparse change detection in
#' high-dimensional linear regression. arXiv preprint, arXiv:2208.06326.
#' @importFrom evd fgev qgev
#' @importFrom RSpectra svds
#' @importFrom MASS mvrnorm
#' @importFrom stats rt rnorm mad predict quantile rbinom rexp rt setNames coef
NULL
