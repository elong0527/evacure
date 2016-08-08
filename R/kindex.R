library(Rcpp)

cppFunction('double varterm_2(NumericMatrix u, NumericVector pi) {
  int nrow = u.nrow();
  double var_k = 0;
  double a = 0;
  for(int i = 0; i < nrow; i++ ){
    for(int j = 0; j < nrow; j++ ){
        if ( i < j ) {
        a += u(j,i) * pi(i) * pi(j)  + u(i,j) * pi(i) * pi(j);
      }
    }
  }
  a = a / nrow / (nrow - 1) * 2;

  for(int i = 0; i < nrow; i++ ){
    for(int j = 0; j < nrow; j++ ){
      for(int k = 0; k < nrow; k++ ){
        if ( i < j & i < k & k != j) {
          var_k +=  (u(j,i) * pi(i) * pi(j)  + u(i,j) * pi(i) * pi(j)  - a) *
                    (u(k,i) * pi(i) * pi(k)  + u(i,k) * pi(i) * pi(k)  - a) ;
        }
      }
    }
  }
  return var_k;
}')

cppFunction('double varterm_3(NumericVector pi) {
  int nrow = pi.length();
  double var_k = 0;
  double a = 0;
  for(int i = 0; i < nrow; i++ ){
    for(int j = 0; j < nrow; j++ ){
      if ( i < j ) {
        a +=  pi(i) * pi(j) ;
      }
    }
  }
  a = a / nrow / (nrow - 1) *2 ;

  for(int i = 0; i < nrow; i++ ){
    for(int j = 0; j < nrow; j++ ){
      for(int k = 0; k < nrow; k++ ){
        if ( i < j & i < k & k != j) {
          var_k +=  ( pi(i) * pi(j)    - a) *
            (pi(i) * pi(k)    - a) ;
        }
      }
    }
  }
  return var_k;
}')

cppFunction('double covterm_2(NumericMatrix u, NumericVector pi) {
  int nrow = u.nrow();
  double var_k = 0;
  double a = 0;
  double b = 0;
  for(int i = 0; i < nrow; i++ ){
    for(int j = 0; j < nrow; j++ ){
      if ( i < j ) {
        a += u(j,i) * pi(i) * pi(j)  + u(i,j) * pi(i) * pi(j);
        b += pi(i) * pi(j) ;
      }
    }
  }
  a = a / nrow / (nrow - 1) * 2  ;
  b = b / nrow / (nrow - 1) * 2 ;

  for(int i = 0; i < nrow; i++ ){
    for(int j = 0; j < nrow; j++ ){
      for(int k = 0; k < nrow; k++ ){
        if ( i < j & i < k & k != j) {

          var_k +=  (u(j,i) * pi(i) * pi(j)  + u(i,j) * pi(i) * pi(j)  - a) *
                    ( pi(i) * pi(k)    - b) ;
        }
      }
    }
  }
  return var_k;
}')

#------------------ Caculate K-index---------------------------------#

.fun <- function(x1,x2, model){  # Kernel Function
  diff <- x1 - x2

  if( model == "PH") weight <- 1/(1+exp(diff))
  if( model == "PO"){ weight <- (1 - exp(diff) + diff * exp(diff)) / (1 - exp(diff))^2
                      weight <- ifelse(is.nan(weight), 0.5, weight)
  }
  if( model == "Normal") weight <- 1 - pnorm(diff/ sqrt(2) )
  (diff < 0) * weight + 0.5 * (diff == 0)
}

#' K-index estimation
#'
#' @export
#' @param risk the risk score for survival part
#' @param pi the uncured probability
#' @param model the type of survival model ("PH", "PO","Normal")
k.ind <- function(risk, pi, model){
  risk <- as.numeric(risk)
  pi   <- as.numeric(pi)
  N <- sum(pi)
  comp.risk.cox <- outer( risk, risk, .fun, model = model)
  diag(comp.risk.cox) <- 0
  k <- sum(comp.risk.cox *   outer(pi, pi, "*")  ) * 2  / ( N^2 - sum(pi^2)) # Based on definition
  k
}

#' K-index estimation (aft model)
#'
#' @export
#' @param risk the risk score for survival part
#' @param pi the uncured probability
#' @param F_est estimated error distribution
k.ind.est <- function(risk, pi, F_est){
  risk <- as.numeric(risk)
  pi   <- as.numeric(pi)
  N <- sum(pi)

  r <- seq(-8,8, length = 2000)
  G <- sapply(r, G_est, F_est = F_est )
  G.foo <- stepfun(r, c(G,0))

  .fun.est <- function(x1, x2, G.foo){
    diff <- x1 - x2
    weight <- G.foo(diff)
    (diff < 0) * weight + 0.5 * (diff == 0) * weight
  }

  # .fun.est <- function(x1, x2, F_est, u){ # Estimated Kernel Function
  #   diff <- x1 - x2
  #   weight <- sapply(diff, G_est, F_est = F_est, u = sort(u) )
  #   (diff < 0) * weight + 0.5 * (diff == 0) * weight
  # }

  comp.risk.cox <- outer( risk, risk, .fun.est, G.foo = G.foo)
  #   comp.risk.cox <- outer( risk, risk, .fun.est, F_est = F_est, u = u)
  diag(comp.risk.cox) <- 0
  k <- sum(comp.risk.cox *   outer(pi, pi, "*")  ) * 2  / ( N^2 - sum(pi^2)) # Based on definition
  k
}

#' Estimate the kernel function G
#'
#' @param F_est CDF of error distribution
#' @param r risk difference
#' @param u the design sequence point
G_est <- function(F_est, r, u = seq(-10,10, length = 2000)){
  a <- sum(  F_est(u + r)  * c( diff( F_est(u) ), 0 ) )
  b <- sum(  F_est(u - r)  * c( diff( F_est(u) ), 0 ) )
  b / (a + b)
}


#' K-index estimation with variance
#'
#' @export
#' @param beta parameters for survival part
#' @param b parameters for cure part
#' @param coef.var estimated variance-covariance matrix for parameters
#' @param X covariates for survival part
#' @param Z covariates for cure part
#' @param model the type of survival model ("PH", "PO","Normal")
k.index <- function(beta, b, coef.var, X, Z, model){
  logit.inv <- function(x) 1 / (1 + exp(-x))
  risk <- as.numeric(X %*% beta)
  pi   <- logit.inv( as.numeric(Z %*% b) )
  coef <- c(beta, b)
  nbeta <- 1:length(beta)
  nb    <- seq(length(beta)+1, length(coef) )

  k_ind <- k.ind(risk, pi, model)


  # Variance Part 1
  var.k1 <- 0
  n <- sum(pi)
  N <- length(pi)
  #     weight <- n * (n-1)  # Based on sampling
  weight <- n^2 - sum(pi^2)           # Based on definition
  weight1 <- N * (N-1)
  B <- 2 * weight / weight1
  A <- B * k_ind
  u <- outer(risk, risk, .fun, model = model)
  diag(u) <- 0

  var.k1_1 <- varterm_2(u, pi) * 4 / weight1^2
  var.k2_1 <- varterm_3(pi)    * 4 / weight1^2
  cov.k1_1 <- covterm_2(u, pi) * 4 / weight1^2
  var.k1 <- 1/B^2 * var.k1_1 + A^2 / B^4 * var.k2_1 - 2 *  A / B^3 * cov.k1_1

  # Variance Part 2
  sd.coef <- sqrt(diag(coef.var)) * 1e-5

  coef.sd    <- coef  + diag(sd.coef)
  coef.sd.ng <- coef  - diag(sd.coef)
  r1 <- apply(coef.sd[nbeta, ], 2, function(xbeta){
    as.numeric(X %*% xbeta )
  })
  r2 <- apply(coef.sd.ng[nbeta, ], 2, function(xbeta){
    as.numeric(X %*% xbeta )
  })
  pi1 <- apply(coef.sd[nb, ], 2, function(xb){
    logit.inv( as.numeric(Z %*% xb ) )
  })
  pi2 <- apply(coef.sd.ng[nb, ], 2, function(xb){
    logit.inv( as.numeric(Z %*% xb ) )
  })


  N.x <- nrow(X)
  D <- ( apply(rbind(r1,pi1),2, function(x) k.ind(x[1:N.x], x[-(1:N.x)], model)) -
           apply(rbind(r2,pi2),2, function(x) k.ind(x[1:N.x], x[-(1:N.x)], model)) ) / 2 / sd.coef
  var.k <- var.k1 + t(D) %*% coef.var %*% D
  se.k <- as.numeric(sqrt(var.k))
  c(k_ind = k_ind, se = se.k)
}

#' Logit transformation
#' @export
logit <- function(x) log(x / (1 - x))

#' Inverse logit transformation
#' @export
logit.inv <- function(x) 1/ (1 + exp(-x))
