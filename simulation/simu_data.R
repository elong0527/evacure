# Simulate Mixture Cure Model

#' Logit transformation
#' @export
logit <- function(x) log(x / (1 - x))

#' Inverse logit transformation
#' @export
logit.inv <- function(x) 1/ (1 + exp(-x))

#' Simulate data from transformation mixture cure models
#'
#' The model is a AFT model (log transformation) with different error
#' distribution. The model option specify the error distriubtion. PH model is a
#' extrme value distribution; PO model is a standard logistic distribution;
#' Normal model is a standard normal distribution.
#'
#' @export
#' @param N sample size
#' @param c.min Uniform censoring minimal bound
#' @param c.max Uniform censoring maximal bound
#' @param model ("PH","PO","Normal") specify error distribution. See description
#' @param .beta parameters for Cox part
#' @param .gamma parameters for logit part (with intercept)
#' @param share wheter cox and logit part share components
#' @param negtive whether - X %*% .beta + error or X %*% .beta + error.true
#' @examples
#' # N <- 200
#' c.min <- 0
#' c.max <- 25
#' model <- "PH"
#' .beta  <- c(- 2, 2)
#' .gamma <- c(1, 0.5, -0.5)
#' nboot <- 3  # Number of bootstrap samples ( Options for smcure1 )
#' var.smcure <- F # Whether caculate bootstrap variance ( Options for smcure1 )
#' simu.cure(N, c.min, c.max, model, .beta, .gamma)
simu.cure <- function(N, c.min, c.max, model, .beta, .gamma, share = T, negative = T){
  library(survival)
  library(smcure)

#   if(share){
#     X    <- cbind(rnorm(N),  rnorm(N) )
#     Z    <- cbind(1, -X)
#   }else{
#     X    <- cbind(rnorm(N),  rnorm(N) )
#     Z    <- cbind(1, rnorm(N), rnorm(N))
#   }

  if(share){
    X    <- cbind(rnorm(N, sd = 2),  rbinom(N, size = 1, prob = 0.5) )
    Z    <- cbind(1, X)
  }else{
    X    <- cbind(rnorm(N),  rbinom(N, size = 1, prob = 0.5) )
    Z    <- cbind(1, rnorm(N), X[,2])
  }

  #   X    <- cbind( runif(N, -2,2),  runif(N, -2, 2) )
  #   Z    <- cbind(1, -X)

  ##  Risk part

  # Transformation
  H.inv <- function(x) exp(x)     # AFT

  # Error term
  error <- function(N, model, G.inv = NULL){
    unif <- runif(N)
    if(model == "PH"){ error <- log(- log(unif)) }
    if(model == "PO"){ error <- logit(unif)}
    if(model == "Normal"){ error <- rnorm(N)}
    if(! model %in% c("PH","PO","Normal") ){ error <- G.inv(unif)}
    error
  }

  error.true <- error(N, model)
  if(negative){
    failure <- H.inv( - X %*% .beta + error.true )
  }else{
    failure <- H.inv( X %*% .beta + error.true )
  }
  censor  <- runif(N, min = c.min, max = c.max)
  # censor  <- rexp(N, rate = c.max)    # Exponential Censoring
  t.all     <- pmin(failure, censor)
  delta.all <- t.all == failure

  # coxph( Surv(failure, rep(1,N)) ~ X ) # Check result

  ## Logit Part
  pi <- logit.inv( Z %*% .gamma )
  cure <- rbinom(N, size = 1,  1 - pi)
#   censor.cure <- runif(N, min = 0, max = 25) # Provides censortime for cured patients

  t      <- ifelse(cure, censor, t.all)
  delta  <- ifelse(cure, 0, delta.all)

  cen.pct.cure <- 1 - sum(delta) / N
  cure.rate <- sum(cure) / N

  info <- c(cen.pct.cure, cure.rate)
  names(info) <- c("cen.rate","cure.rate")

  surv_res <- cbind(t,delta,cure,failure)
  colnames(surv_res) <- c("t","delta","cure","failure")

  data <- list(surv = surv_res,
               X = X,
               Z = Z,
               beta = .beta,
               gamma = .gamma,
               info = info
  )

  # Check result
  # smcure( Surv(t, delta) ~  I(- X.1) + I(- X.2) ,
  #         cureform = ~  Z.1 + Z.2, model = "aft",
  #         data = data, Var = var.smcure, nboot = nboot)

}

