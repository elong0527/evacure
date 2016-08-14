

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


  comp.risk.cox <- outer( risk, risk, .fun.est, G.foo = G.foo)
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

#' Logit transformation
#' @export
logit <- function(x) log(x / (1 - x))

#' Inverse logit transformation
#' @export
logit.inv <- function(x) 1/ (1 + exp(-x))
