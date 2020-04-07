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
  library(mvtnorm)
  # print(nboot)
  simu <- simu.cure(N, c.min, c.max, model, .beta, .gamma, share = share)
  data <- data.frame(simu$surv, X = simu$X, Z = simu$Z[,-1])

  fit <- smcure1( formula = Surv(t, delta) ~  X.1 + X.2 ,
                  cureform = ~  Z.1 + Z.2, model = fit_method,
                  data = data, Var = var.smcure, nboot = nboot, init = init, cutpoint = cutpoint)

  surv_model    <- model.frame(fit$call$formula, data = data)

  # Direct estimate spe and sen
  X  <- as.matrix(surv_model[,-1])
  beta <- fit$beta
  Z  <- as.matrix(cbind(Intercept = 1, model.frame(fit$call$cureform, data = data )))
  b <- fit$b

  post_fit <- eva_cure_direct(time = data$t, data$delta, X, beta, Z, b, model, cutpoint, n_post = 500)
  direct_est <- cbind(t(post_fit$metric), C = NA, sen = t(post_fit$senspe[,1]), spe = t(post_fit$senspe[,2]))


  if(var.smcure){
    b_boot    <- lapply(fit$eva.boot, function(x) x$para[-(1:ncol(X))])
    beta_boot <- lapply(fit$eva.boot, function(x) x$para[(1:ncol(X))])

    b_imp_mean <- apply(do.call(rbind, b_boot), 2, mean)
    b_imp_var  <- var((do.call(rbind, b_boot)))

    beta_imp_mean <- apply(do.call(rbind, beta_boot), 2, mean)
    beta_imp_var  <- var((do.call(rbind, beta_boot)))

    # Direct use Bootstrap Sample
    direct_est_boot1 <- list()
    for(i_boot in 1:nboot){
      .boot1 <- eva_cure_direct(time = data$t, data$delta, X, beta = beta_boot[[i_boot]], Z, b = b_boot[[i_boot]], model, cutpoint, n_post = 500)
      .direct_est1 <- cbind(t(post_fit$metric), C = NA, sen = t(.boot1$senspe[,1]), spe = t(.boot1$senspe[,2]))
      direct_est_boot1[[i_boot]] <- c(b = b_boot[[i_boot]], beta = beta_boot[[i_boot]], .direct_est1)
    }
    direct_est_boot1 <- do.call(rbind, direct_est_boot1)

    # Use sample from normal distribution
    direct_est_boot2 <- list()
    for(i_boot in 1:nboot){

      b_imp <- rmvnorm(1, b_imp_mean, b_imp_var)
      beta_imp <- rmvnorm(1, beta_imp_mean, beta_imp_var)

      .boot2 <- eva_cure_direct(time = data$t, data$delta, X, beta = beta_boot[[i_boot]], Z, b = b_boot[[i_boot]], model, cutpoint, n_post = 500)
      .direct_est2 <- cbind(t(post_fit$metric), C = NA, sen = t(.boot2$senspe[,1]), spe = t(.boot2$senspe[,2]))
      direct_est_boot2[[i_boot]] <- c(b = b_boot[[i_boot]], beta = beta_boot[[i_boot]], .direct_est2)
    }
    direct_est_boot2 <- do.call(rbind, direct_est_boot2)


  }

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
  est  = c(fit$b, fit$beta, fit$eva),
  direct = c(b = simu$gamma, beta = simu$beta, direct_est)
  )

  if(var.smcure == T){
    b_sd0    <- apply(fit$b_boot, 2, function(x) c( mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975))) )
    beta_sd0 <- apply(fit$beta_boot, 2, function(x) c( mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975))) )

    boot1_res <- apply(direct_est_boot1, 2, function(x) c( boot1.mean = mean(x), boot1.sd = sd(x), boot1 = quantile(x, c(0.025, 0.975), na.rm = TRUE)) )
    boot2_res <- apply(direct_est_boot2, 2, function(x) c( boot2.mean = mean(x), boot2.sd = sd(x), boot2 = quantile(x, c(0.025, 0.975), na.rm = TRUE)) )

    res <- rbind(res,
                 cbind(b_sd0, beta_sd0, fit$eva_sd),
                 boot1_res,
                 boot2_res)

    coef_var <- var(cbind(fit$beta_boot, fit$b_boot))
    res <- res
    #   res <- round(res, 2)
  }


  list(res = res, info = simu$info)
}

