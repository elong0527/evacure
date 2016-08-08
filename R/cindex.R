#' C index for mixture cure models
#'
#' Nonparametrics estimators for C index
#' @export
#' @param risk the risk score for survival part
#' @param pi the uncured probability
#' @param time observation time
#' @param delta censoring status
#' @param tau boundry of time
#' @examples
#' fit <- coxph(Surv(time, status) ~  age + sex , lung, x = T)
#' risk <- fit$x %*% fit$coefficients
#' survConcordance(Surv(time, status) ~ risk, lung)
#' pi <- rep(1,fit$n)
#' c.index(risk, pi, lung$time, lung$status ==2)
#
#
#' fit <- coxph(Surv(stop, event) ~  rx + size , bladder, x = T)
#' risk <- fit$x %*% fit$coefficients
#' survConcordance(Surv(stop, event) ~ risk, bladder)
#' pi <- rep(1,fit$n)
#' c.index(risk, pi, bladder$stop, bladder$event)
cindex <- function(risk, pi, time, delta, tau = NULL){
  risk <- as.numeric(risk)
  pi   <- as.numeric(pi)
  n <- length(risk)
  t1 <- outer(time, time, "<"  )
  t2 <- rep(1, n)
  if(! is.null(tau) )t2 <- time < tau  # Boundry of time
  r1 <- outer(risk, risk, function(x,y) (x>y) + 0.5 * (x==y))
  p1 <- outer(pi, pi ,"*")
  num <- delta * t1 * t2 * r1 * p1
  diag(num) <- 0
  dev <- delta* t1 * t2 * p1
  diag(dev) <- 0
  sum(num) / sum(dev)
}





