getIntervals <- function(n, p, no_intervals, delta=0.1){
  intervals <- matrix(sample.int(n, no_intervals*100, replace=TRUE), ncol = 2)
  intervals <- t(apply(intervals, 1, sort))
  intervals <- intervals[intervals[,2] - intervals[,1] >= (1 + delta)*p, ]
  if (nrow(intervals) > no_intervals) intervals <- intervals[1:no_intervals, ]
}

#' Generate the threshold for Algorithm 3 in (Gao and Wang 2022)
#' @param n the sample size, i.e., the number of rows in the design
#' @param p the regression coefficient dimension 
#' @param burnIn the burnIn parameter to be passed the cpreg function
#' @param alpha the critical type-I error threshold, should in practice adjusted for the multiple testing, i.e., in this case the numer of intervals
#' @param sigma the standard deviation of the noise in the regressions  
#' @param b the standard deviation of each entry of the randomly generated regression coeffients, which follows N(0,b^2)
#' @return a nonnegative scalar as the threshold for not_cpreg
#' @description Generate the threshold for combining charcoal with the narrowest-over-threshold algorithm, which is also used as the critical threshold in the testing refining stage in not_cpreg.  It works by generating permSize Monte Carlo repetetions of running cpreg on the null model without any change in regression coefficients 
#' @export
getThreshold <- function(n, p, burnIn=0, alpha=0.05, permSize=1000, sigma=1,
                         b=1) {
  cpreg_stats <- rep(-Inf, permSize)
  for (i in 1:permSize) {
    X <- matrix(rnorm(n*p), nrow=n)
    beta <- rnorm(p) * max(b,1)
    Y <- X %*% beta + rnorm(n) * sigma
    tmp <- cpreg(X, Y, sigma=sigma, burnIn=burnIn)
    cpreg_stats[i] <- tmp$stats[tmp$cp]
  }

  if (alpha * permSize >= 1){
    zeta <- quantile(cpreg_stats, 1 - alpha)
  } else {
    gev.fit <- evd::fgev(cpreg_stats)$param
    zeta <- evd::qgev(1 - alpha, loc=gev.fit[1], scale=gev.fit[2],
                      shape=gev.fit[3])
  }
  return(list(thresh=zeta, stats=cpreg_stats))
}

#' Soft thresholding a vector
#' @param x a vector of real numbers
#' @param lambda soft thresholding value
#' @return a vector of the same length
#' @description entries of v are moved towards 0 by the amount lambda until they hit 0.
vector.soft.thresh <- function(x, lambda){
  sign(x)*pmax(0,(abs(x)-lambda))
}


#' Compute matrix power of a symmetric matrix
#' @param A a square matrix
#' @param power power exponent
#' @param pseudoinverse whether to use pseudoinverse if power is negative
matrix.power <- function(A, power, pseudoinverse=TRUE){
  if (nrow(A)!=ncol(A)) stop('A need to be a square matrix.')
  symm <- isSymmetric(A)
  tmp <- eigen(A, symmetric=symm)
  evecs <- tmp$vectors
  evals <- tmp$values
  evals[abs(evals) < 1e-12] <- 0
  
  if (sum(evals==0) > 0 && power < 0){
    if (!pseudoinverse){
      stop()
    } else {
      power <- -power
      evals[evals!=0] <- 1/evals[evals!=0]
    }
  }
  
  if (symm && min(Re(evals)) >=0) {
    return(evecs %*% diag(evals^power, nrow=nrow(A)) %*% t(evecs))
  } else {
    return(evecs %*% diag(as.complex(evals)^power, nrow=nrow(A)) %*% t(evecs))
  }
}

#' Print percentage
#' @param ind a vector of for loop interator
#' @param tot a vector of for loop lengths
#' @return on screen output of percentage
printPercentage <- function (ind, tot){
  ind <- as.vector(ind); tot <- as.vector(tot)
  if ((length(tot) > 1) & (length(ind) == 1)) {ind <- match(ind, tot); tot <- length(tot)}
  len <- length(ind)
  contrib <- rep(1,len)
  if (len > 1) {
    for (i in (len-1):1) contrib[i] <- contrib[i+1] * tot[i+1]
  }
  grand_tot <- contrib[1] * tot[1]
  count <- (sum(contrib * (ind - 1)) + 1)
  out <- ""
  if (sum(ind-1)>0) out <- paste0(rep("\b", nchar(round((count-1)/grand_tot * 100))+1), collapse = "")
  out <- paste0(out, round(count/grand_tot*100), "%")
  if (identical(ind, tot)) out <- paste0(out, '\n')
  cat(out)
  return(NULL)
}

