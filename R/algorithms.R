library(glmnet)



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



#' Computing the orthonormal matrix spanning the orthogonal complement of the column span of X
#' @param X an n x p matrix with n > p
#' @return matrix A of size n x (n-p) with orthonormal columns
#' @description it seems that QR decomposition in R is slow when the dimension is large e.g. n and p > 5000; when working with large matrices, we recommend use either the python or matlab version of this package.
#' @export
orthogonalProjection <- function(X){
  n <- dim(X)[1]; p <- dim(X)[2]; m <- n - p
  tmp <- matrix(rnorm(n*m), n, m)
  tmp <- tmp - X%*%(solve(t(X)%*%X)%*%(t(X)%*%tmp))
  return(qr.Q(qr(tmp)))
}

#' Computing the noise variance in a high-dimensional linear model
#' @param W design matrix for high-dimensional linear model
#' @param z response in a high-dimensional linear model
#' @return a nonnegative scalar of estimated noise standard deviation
#' @description Assume data are generated from the model z = W beta + e, where e has independent components N(0, sigma^2), we estimate sigma^2 by fitting a cross-validated Lasso estimator of beta and then estimate sigma from the residual sum of squares.
#' @export
noise_sd <- function(W, z){
  fit.cv <- glmnet::cv.glmnet(W,z,intercept=F,standardize=F)
  betahat <- coef(fit.cv, s='lambda.min')
  Xbetahat <- predict(fit.cv, newx=W, s='lambda.min')
  sqrt(sum((z - Xbetahat)^2) / (nrow(W) - sum(betahat!=0)))
}

#' Computing the noise variance in a high-dimensional linear model
#' @param W design matrix for high-dimensional linear model
#' @param z response in a high-dimensional linear model
#' @return a nonnegative scalar of estimated noise standard deviation
#' @description Assume data are generated from the model z = W beta + e, where e has independent components N(0, sigma^2), we estimate sigma^2 by using the method in Dicker (2014).
#' @references Dicker, L. H. (2014). Variance estimation in high-dimensional linear models. Biometrika, 101(2), 269-284.
#' @export
dicker_noise_sd <- function(W, z)
{
  m <- nrow(W); p <- ncol(W)
  gram.W.norm <- t(W) %*% W/m
  m1.hat <- sum(diag(gram.W.norm))/p
  m2.hat <- sum(gram.W.norm^2)/p - p * m1.hat^2 / m
  sigma.tilde.square <- (1 + p*m1.hat^2/(m+1)/m2.hat)*sum(z^2)/m - m1.hat*sum((t(W) %*% z)^2)/m/(m+1)/m2.hat
  if (sigma.tilde.square > 0){
    return(sqrt(sigma.tilde.square))
  }
  else return(noise_sd(W,z))
  return(sqrt(sigma.tilde.square))
}


#' Main function implementing the charcoal algorithm in Gao and Wang (2022)
#' @param X design matrix of the linear regression model
#' @param y response vector of the linear regression model
#' @param sparse whether to test a sparse difference in the regression coefficients of the two models, default is TRUE
#' @param sigma noise standard deviation of the regression models; if unknown, set to NULL and will be estimated.
#' @param aggregate_method specifies which variant of the algorithms to be applied. Set it to 1 for Algorithm 1 (default), and to 2 for Algorithm 2.
#' @return a list containing the test statistics and changepoint estimate.
#' @description From the model y_t = x_t beta_1 + eps_t for 1<=t<=z and y_t = x_t beta_2 + eps_t for z+1<=t<=n, where x_t are p-dimensional design vectors for 1<=t<=n, we estimate the changepoint z where the change of regression coefficient takes place. The flag sparse is set to be TRUE if beta1 - beta2 has at most sqrt(p) so that appropriate testing statistics are formed. The test is performed using the charcoal algorithm from Gao and Wang (2022). 
#' @references Gao, F. and Wang, T. (2022) Sparse change detection in high-dimensional regression. Work-in-progress.
#' @export
#' @examples
#' # problem parameters
#' n <- 240 # total sample size
#' z <- 150  # changepoint location
#' p <- 160 # dimension of covariates
#' k <- 10 # sparsity of difference of the two regression coefficients
#' rho <- 1 # difference in l_2 norm of the two regression coefficients
#' # generate design matrices
#' X1 <- matrix(rnorm(z * p), z, p)
#' X2 <- matrix(rnorm((n-z) * p), (n-z), p)
#' # generate regression coefficients
#' beta1 <- rnorm(p)
#' theta <- c(rnorm(k), rep(0, p-k))
#' theta <- theta / sqrt(sum(theta^2)) * rho
#' beta2 <- beta1 + theta
#' # generate response vectors
#' y1 <- X1 %*% beta1 + rnorm(n1)
#' y2 <- X2 %*% beta2 + rnorm(n2)
#' # combine the samples before and after changepoint to form the response and covariates
#' X <- rbind(X1,X2)
#' Y <- rbind(y1,y2)
#' # test for difference in beta1 and beta2
#' result <- cpreg(X,Y)
#' # or  result <- cpreg(X,Y,lambda=NA,sigma=NULL,aggregate_method=2)
cpreg <- function(X, Y, lambda=NA, sigma=NULL, aggregate_method=1){
  n <- dim(X)[1]; p <- dim(X)[2]; m <- n - p
  if (is.na(lambda)) lambda <- 2*sqrt(log(p))
  if (is.null(sigma)) sigma_recalc <- TRUE else sigma_recalc <- FALSE
  
  A <- t(orthogonalProjection(X))
  Z <- A %*% Y
  
  W <- matrix(0, m, p)
  Q <- matrix(0, p, n-1)
  Q2 <- matrix(0, p, n-1)
  for (t in 1:(n-1)) {
    W <- W + 2 * outer(A[,t], X[t,])
    if (sigma_recalc)  sigma = dicker_noise_sd(W,Z)
    W = W / sigma; Z = Z/sigma
    col_scale <- sqrt(colSums(W^2))
    W_tilde = sweep(W, 2, col_scale, FUN="/")
    Q[, t] <- t(W_tilde) %*% Z
    if (aggregate_method==3) Q2[, t] <- matrix.power(t(W)%*%W, -1/2) %*% t(W) %*% Z
    printPercentage(t, n-1)
  }
  stats <- sqrt(colSums(vector.soft.thresh(Q, lambda)^2))
  if (aggregate_method==1){
    stats <- sqrt(colSums(vector.soft.thresh(Q, lambda)^2))
  } else if (aggregate_method==2){
    Q.thresh <- vector.soft.thresh(Q, lambda)
    v <- svd(Q.thresh)$v[,1]
    stats <- abs(v)
  } else if (aggregate_method==2.5){
    Q.thresh <- vector.soft.thresh(Q, lambda)
    u <- svd(Q.thresh)$u[,1]
    stats <- as.vector(abs(t(u)%*%Q))
  } else if (aggregate_method==0){
    stats <- 0
  } else if (aggregate_method==3){
    stats <- sqrt(colSums(vector.soft.thresh(Q2, lambda)^2))
  }
  cp <- which.max(stats)

  return(list(cp=cp, stats=stats, Q=Q))
}
