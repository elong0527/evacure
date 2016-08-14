#' Simulation for coxph cure model
#'
#' Report Sensitivity/Specificity AUC C-index K-index for
#' coxph mixture cure model
#'
#' @export
#' @param N sample size
#' @param c.min Uniform censoring minimal bound
#' @param c.max Uniform censoring maximal bound
#' @param model ("PH","PO","Normal") specify error distribution. See description
#' @param .beta parameters for Cox part
#' @param .gamma parameters for logit part (with intercept)
#' @param share wheter cox and logit part share components
#' @param fit_method the method to fit the cured model ("ph" or "aft")
#' @examples
# N <- 200
# c.min <- 0
# c.max <- 25
# model <- "PH"
# .beta  <- c(- 2, 2)
# .gamma <- c(1, 0.5, -0.5)
# nboot <- 3  # Number of bootstrap samples ( Options for smcure1 )
# var.smcure <- T # Whether caculate bootstrap variance ( Options for smcure1 )
# share = T
# fit_method = "ph"
# init = NULL
# cutpoint = c(0.1,0.25,0.5,0.75,0.9)
simu.coxph <- function(N, c.min, c.max, model, .beta, .gamma, share = T, var.smcure = F, fit_method = "ph",
                         init = NULL, cutpoint = c(0.1,0.25,0.5,0.75,0.9), nboot = 50 ){
  options(warn = -1)
  library(ROCR)
  # print(nboot)
  simu <- simu.cure(N, c.min, c.max, model, .beta, .gamma, share = share)
  data <- data.frame(simu$surv, X = simu$X, Z = simu$Z[,-1])

  fit <- smcure1( formula = Surv(t, delta) ~  X.1 + X.2 ,
                  cureform = ~  Z.1 + Z.2, model = fit_method,
                  data = data, Var = var.smcure, nboot = nboot, init = init, cutpoint = cutpoint)

  # True metrics
  i.uncure <- data$cure == 0
  tr.odds0 <- simu$Z %*% simu$gamma
  tr.risk0 <- simu$X %*% simu$beta
  tr.risk <-  tr.risk0[i.uncure]
  tr.odds <- tr.odds0[i.uncure]
  tr.pi   <- rep(1, sum(i.uncure))
  tr.t    <- data$t[i.uncure]
  tr.delta <- data$delta[i.uncure]
  tr.pi0   <- logit.inv(tr.odds0)
  tr.surv0  <- exp(- data$t * exp( tr.risk0 ) )
  tr.w0  <- w.cure(tr.pi0, data$delta, tr.surv0)

  k0 <- k.ind(risk = tr.risk, pi = tr.pi, model = model)
  c0 <- cindex(risk = tr.risk, pi = tr.pi, time = tr.t, delta = tr.delta )
  a0 <- performance( prediction(tr.odds0, 1 - data$cure ), measure = "auc")@y.values[[1]]
  a01<- auc.cure(tr.odds0, tr.risk0, tr.surv0, data$delta)

  sen <- performance( prediction(tr.odds0, 1 - data$cure ),  "sens" )
  spe <- performance( prediction(tr.odds0, 1 - data$cure ),  "spec" )
  res0 <- cbind(sen = sen@y.values[[1]], spe = spe@y.values[[1]], cut = spe@x.values[[1]])
  res01<- cutSenspe(res0, cutpoint )

  res <- rbind(
  true = c(b = simu$gamma, beta = simu$beta, AUC = a01, K = k0, C = c0, sen = res01[,1], spe = res01[,2]),
  est  = c(fit$b, fit$beta, fit$eva)
  )

  if(var.smcure == T){
    res <- rbind(res,   cbind(fit$b_sd, fit$beta_sd, fit$eva_sd) )

    coef_var <- var(cbind(fit$beta_boot, fit$b_boot))
    res <- res
    #   res <- round(res, 2)
  }

  list(res = res, info = simu$info)
}
