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
#' @description Generate the orthogonal complement of the column span of matrix X 
#' @export
orthogonalProjection <- function(X){
  n <- dim(X)[1]; p <- dim(X)[2]
  qr.Q(qr(X), complete=TRUE)[, (p+1):n]
}

#' Computing the noise variance in a high-dimensional linear model
#' @param W design matrix for high-dimensional linear model
#' @param z response in a high-dimensional linear model
#' @return a nonnegative scalar of estimated noise standard deviation
#' @description Assume data are generated from the model z = W b + e, where e has independent components N(0, s^2), we estimate s^2 by fitting a cross-validated Lasso estimator of b and then estimate sigma from the residual sum of squares.
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
#' @description Assume data are generated from the model z = W b + e, where e has independent components N(0, s^2), we estimate s^2 by using the method in Dicker (2014).  If the method fails, a fallback with nosie_sd will be used.
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
  } else return(noise_sd(W,z))
}





#' Main function implementing the charcoal algorithms for single changepoints
#' in Gao and Wang (2022), including Algorithms 1, 2 and 3
#' @param X design matrix of the linear regression model
#' @param Y response vector of the linear regression model
#' @param lambda threshold for Algorithms 1 and 2, set it to NULL and the recommended value will be used 
#' @param sigma noise standard deviation of the regression models; if unknown, set to NULL.
#' @param aggregate_method specifies which variant of the algorithms to be applied. Set it to 'compsket' for Algorithm 1, to 'proj' for Algorithm 2, and to 'lasso_bic' for Algorithm 3.
#' @param burnIn specifies the fraction at both ends of the interval to be discarded as possbile changes, to handle common boundary effects
#' @return a list containing the test statistics and changepoint estimate.
#' @description From the model y_t = x_t beta_1 + eps_t for 1<=t<=z and y_t = x_t beta_2 + eps_t for z+1<=t<=n, where x_t are p-dimensional design vectors for 1<=t<=n, we estimate the changepoint z where the change of regression coefficient takes place. The localization is performed using the charcoal algorithms from Gao and Wang (2022). 
#' @references Gao, F. and Wang, T. (2022) Sparse change detection in high-dimensional linear regression.
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
#' # combine the samples before and after changepoint to form the
#' # response and covariates
#' X <- rbind(X1,X2)
#' Y <- rbind(y1,y2)
#' # test for difference in beta1 and beta2
#' result <- cpreg(X,Y)
#' # or  result <- cpreg(X,Y) 
cpreg <- function(X, Y, lambda=NA, sigma=NULL, burnIn=0, aggregate_method='proj')
{
  n <- dim(X)[1]; p <- dim(X)[2]; m <- n - p
  if (p>=n) return(list(stats=NULL,cp=NULL))
  if (!is.null(sigma)) {X <- X/sigma; Y<- Y/sigma}

  if (is.na(lambda)) lambda <- 0.5*log(p) # default recommended setting

  A <- t(orthogonalProjection(X))
  Z <- A %*% Y

  W <- matrix(0, m, p)
  Q <- matrix(0, p, n-1)
  stats <- rep(-Inf, n-1)
  coeffs <- matrix(0, p, n-1)

  if (grepl('lasso',aggregate_method)) {
    for (t in 1:(n-1)) {
      W <-  W + 2 * outer(A[,t], X[t,])
      if (t >= 5 && t <= n-5){
        tmp <- glmnet::cv.glmnet(W, Z,nfolds=5, intercept=FALSE,
                               standardize=TRUE, parallel = FALSE)
        b <- as.vector(coef(tmp, s=tmp$lambda.min))[-1]
        coeffs[, t] <- b
        if (grepl('bic',aggregate_method)) {
          stats[t] <- -(sum((Z - W %*% b)^2) + log(m) * sum(b!=0))
        } else stats[t] <- -(sum((Z - W %*% b)^2))
      }
    } } else {
      col_scales <- matrix(0, p, n-1)
      for (t in 1:(n-1)) {
        W <- W + outer(A[,t], X[t,])
        col_scales[,t] <- sqrt(colSums(W^2))
      }
      Q <- X[-n, ] * colSums(A[, -n] * as.vector(Z))
      Q <- t(apply(Q, 2, cumsum)) / col_scales
    }

  idx_begin <- (round(burnIn*n)+1)
  idx_end <- (-round(burnIn*n)+n-1)

  lambda <- lambda * mad(Q)
  if (aggregate_method=='compsket'){
    stats <- sqrt(colSums(vector.soft.thresh(Q[, ], lambda)^2))
  } else if (aggregate_method=='svd'){
    Q.thresh <- vector.soft.thresh(Q[, ], lambda)
    v <- svd(Q.thresh)$v[,1]
    stats <- abs(v)
  } else if (aggregate_method=='proj'){
    Q.thresh <- vector.soft.thresh(Q[, ], lambda)
    # u <- svd(Q.thresh)$u[,1]
    u <- RSpectra::svds(Q.thresh, 1)$u
    stats <- as.vector(abs(t(u)%*%Q[, ]))
  } else if (grepl('lasso', aggregate_method)) {
    attr(stats, 'coeff') <- coeffs[,]
  } else {
    stats <- NA
  }

  cp <- which.max(stats[(idx_begin:idx_end)]) +idx_begin - 1
  return(list(cp=cp, stats=stats, Q=Q, besidesBurned=c(idx_begin,
  idx_end)))
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



#' Main function implementing the charcoal algorithms for single changepoints
#' in Gao and Wang (2022), including Algorithms 1, 2 and 3
#' @param X design matrix of the linear regression model
#' @param Y response vector of the linear regression model
#' @param sigma noise standard deviation of the regression models; if unknown, set to NULL.
#' @param cpreg_method specifies which variant of the algorithms to be applied in the refinement stage. Set it to 'compsket' for Algorithm 1, to 'proj' for Algorithm 2, and to 'lasso_bic' for Algorithm 3. The default is 'lasso_bic'.
#' @param burnIn the burnIn parameter to be passed to cpreg, which specifies the fraction at both ends of the interval to be discarded as possbile changes, to handle common boundary effects.
#' @return a list containing the changepoint estimates, the estimates before testing refinement, the estimates after the midpoint refinement and the estimates after the final refinement.
#' @description From the model y_t = x_t beta_1 + eps_t for 1<=t<=z and y_t = x_t beta_2 + eps_t for z+1<=t<=n, where x_t are p-dimensional design vectors for 1<=t<=n, we estimate the changepoint z where the change of regression coefficient takes place. 
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
#' # combine the samples before and after changepoint to form the
#' # response and covariates
#' X <- rbind(X1,X2)
#' Y <- rbind(y1,y2)
#' # test for difference in beta1 and beta2
#' result <- cpreg(X,Y)
#' # or  result <- cpreg(X,Y) 
not_cpreg <- function(X, Y, zeta=NULL,
                      no_intervals = floor(nrow(X)/5), burnIn=0, sigma=NULL,
                      verbose=FALSE, cpreg_method='lasso_bic', seed=NULL) {
  # zeta --- threshold for narrowest-over-threshold

  n <- nrow(X); p <- ncol(X)
  if (is.null(sigma)) {
    ret <- cpreg(X, Y)
    sigma <- mad(ret$Q)
  }
  if (is.null(zeta)) zeta <- getThreshold(n, p, burnIn, alpha=0.05/n,
                                          1000, sigma = sigma, verbose=verbose)
  if (!is.null(seed)) set.seed(seed)

  # get intervals, left open, right closed intervals
  intervals <- getIntervals(n, p, no_intervals, 0.1)

  not <- function(start, end){
    if (verbose) print(paste0('Searching interval: (', start, ', ', end, ']'))
    if (end - start <= p) return(NULL)

    M_se  <- intervals[intervals[,1] >= start & intervals[,2] <= end, ]
    M_se <- rbind(M_se, c(start, end))
    results <- data.frame(left=M_se[,1], right=M_se[,2])
    results$lengths <- results$right - results$left
    results$cps <- results$stats <- 0

    for (i in 1:nrow(results)){
      sm <- results$left[i]; em <- results$right[i]
      tmp <- cpreg(X[(sm+1):em,], Y[(sm+1):em], sigma=sigma, burnIn=burnIn)
      results$cps[i] <- tmp$cp + sm
      results$stats[i] <- tmp$stats[tmp$cp]
      results$sigma[i] <- mad(tmp$Q)
    }

    if (sum(results$stats >= zeta) == 0) return(NULL)
    results <- results[results$stats >= zeta, ]
    m_star <- which.min(results$lengths)
    cp <- results$cps[m_star]
    if (verbose)  print(c('discovered cp at', cp))

    ret <- setNames(c(cp, results$stats[m_star]),
                    c('cps','stats'))
    return(rbind(not(start, cp), ret, not(cp, end)))
  }

  # Stage 1: initial changepoint searching via NOT
  ret_initial <- not(0, n)
  if (is.null(ret_initial)) return(NULL)

  ret_initial <- as.data.frame(ret_initial, row.names=1:nrow(ret_initial))
  if (verbose) print(ret_initial)

  # Stage 2: perform tests on each identified changepoint
  ret_tested <- ret_initial
  for (i in order(ret_initial$stats)){
    cp <- ret_initial$cps[i]
    before.idx <- max(ret_tested$cps[ret_tested$cps < cp], 0)
    after.idx <- min(ret_tested$cps[ret_tested$cps > cp], n)
    if (verbose) print(paste0('testing cp at ', cp, ' in (',
                              before.idx, ', ', after.idx, ']'))

    too_small <- after.idx - before.idx < p + 10
    too_close <- min(cp - before.idx, after.idx - cp) < burnIn * n
    if (!too_small && !too_close){
      tmp <- cpreg(X[(before.idx+1):after.idx, ],
                   Y[(before.idx+1):after.idx], burnIn = burnIn)
      test_stat <- tmp$stats[tmp$cp] / sigma
      test_fail <- test_stat < zeta
    } else {
      test_fail <- FALSE
    }
    if (too_small || too_close || test_fail) {
      ret_tested <- ret_tested[-match(cp, ret_tested$cps), ]
    } else {
      ret_tested[match(cp, ret_tested$cps), 'stats'] <- test_stat
    }
    if (verbose) {
      if (too_small) {
        print('current interval too small, merged')
      } else if (too_close) {
        print('current cp close to boundary, merged')
      } else if (test_fail){
        print('test stat below threshold, removed')
      } else {
        print(paste0('keeping cp at ', cp))
      }
    }
  }
  if (verbose) print(ret_tested)

  # Stage 3: refine the change point estimate via midpoint refinement
  ret_refined <- ret_tested
  for (i in order(ret_tested$stats, decreasing=TRUE)) {
    cp <- ret_refined$cps[i]
    before.idx <- max(ret_refined$cps[ret_refined$cps < cp], 0)
    after.idx <- min(ret_refined$cps[ret_refined$cps > cp], n)
    before.mid <- floor((cp + before.idx) / 2)
    after.mid <- floor((cp + after.idx) / 2)

    if (verbose) print(paste0('midpoint refinement interval (',
                             before.mid, ', ', after.mid, ']'))

    if (after.mid - before.mid >= (1 + burnIn) * p) {
      tmp <- cpreg(X[(before.mid+1):after.mid, ],
                   Y[(before.mid+1):after.mid],
                   burnIn=0, sigma=sigma, aggregate_method=cpreg_method)
      ret_refined$cps[i] <- tmp$cp + before.mid
      ret_refined$stats[i] <- tmp$stats[tmp$cp]
      if (verbose) print(paste0('refining cp at ', cp, ' to new cp at ',
                               ret_refined$cps[i]))
    } else {
      ret_refined$stats[i] <- 0
      if (verbose) print('not big enough interval, skip first refinement')
    }
  }
  if (verbose) print(ret_refined)

  # Stage 4: second refinement
  ret_final <- ret_refined

  for (i in order(ret_tested$stats, decreasing=TRUE)) {
    cp <- ret_final$cps[i]
    burnInLen <- min(10, floor(burnIn*n))
    before.idx <- max(ret_final$cps[ret_final$cps < cp] + burnInLen, 0)
    after.idx <- min(ret_final$cps[ret_final$cps > cp] - burnInLen, n)
    if (after.idx - before.idx >= (1+burnIn) * p) {
      tmp <- cpreg(X[(before.idx+1):after.idx,],
                   Y[(before.idx+1):after.idx],
                   burnIn=burnIn, sigma=sigma, aggregate_method = cpreg_method)
      ret_final$cps[i] <- tmp$cp + before.idx
      ret_final$stats[i] <- tmp$stats[tmp$cp]
      if (verbose) print(paste('refining cp at', cp, 'to', ret_final$cps[i]))
    } else {
      ret_final$stats[i] <- 0
      if (verbose) print('not big enough interval, skip second refinement')
    }
  }
  if (verbose) print(ret_final)

  # return(list(S=S,Stats=Stats,Starts=Starts,Ends=Ends))
  return(list(cps=ret_final$cps, initial=ret_initial, tested=ret_tested,
              refined=ret_refined, final=ret_final))
}
