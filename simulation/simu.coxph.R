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
                         init = NULL, cutpoint = c(0.1,0.25,0.5,0.75,0.9), nboot = 100, n_post = 100 ){
  options(warn = -1)
  library(ROCR)
  library(mvtnorm)
  # print(nboot)
  simu <- simu.cure(N, c.min, c.max, model, .beta, .gamma, share = share)
  simu_test <- simu.cure(N, c.min, c.max, model, .beta, .gamma, share = share)

  data <- data.frame(simu$surv, X = simu$X, Z = simu$Z[,-1])

  fit <- smcure1( formula = Surv(t, delta) ~  X.1 + X.2 ,
                  cureform = ~  Z.1 + Z.2, model = fit_method,
                  data = data, Var = var.smcure, nboot = nboot, n_post = n_post, init = init, cutpoint = cutpoint)

  fit_em_like <- smcure1( formula = Surv(t, delta) ~  X.1 + X.2 ,
                          cureform = ~  Z.1 + Z.2, model = fit_method,
                          data = data, Var = var.smcure, nboot = nboot, n_post = n_post, cutpoint = cutpoint,
                          init = list(b = fit$b, beta = fit$beta, w = fit$w), est_type = "EM-like")

  fit_b <- smcure1( formula = Surv(t, delta) ~  X.1 + X.2 ,
                    cureform = ~  Z.1 + Z.2, model = fit_method,
                    data = data, Var = var.smcure, nboot = nboot, n_post = n_post, init = init, cutpoint = cutpoint, posterior = FALSE)


  n_x <- ncol(fit$X)
  post_fit <- fit$eva_direct
  direct_est <- cbind(t(post_fit$metric), C = NA, sen = t(post_fit$senspe[,1]), spe = t(post_fit$senspe[,2]))

  n_x <- ncol(fit_b$X)
  post_fit <- fit_b$eva_direct
  direct_est_b <- cbind(t(post_fit$metric), C = NA, sen = t(post_fit$senspe[,1]), spe = t(post_fit$senspe[,2]))

  if(var.smcure){
    b_boot    <- lapply(fit$eva.boot, function(x) x$para[-(1:n_x)])
    beta_boot <- lapply(fit$eva.boot, function(x) x$para[(1:n_x)])

    # Direct use Bootstrap Sample

    direct_est_boot1 <- list()
    for(i_boot in 1:nboot){
      .boot1 <- fit$eva_boot_direct[[i_boot]]
      .direct_est1 <- cbind(t(.boot1$metric), C = NA, sen = t(.boot1$senspe[,1]), spe = t(.boot1$senspe[,2]))
      direct_est_boot1[[i_boot]] <- c(b = b_boot[[i_boot]], beta = beta_boot[[i_boot]], .direct_est1)
    }
    direct_est_boot1 <- do.call(rbind, direct_est_boot1)

    b_boot    <- lapply(fit_b$eva.boot, function(x) x$para[-(1:n_x)])
    beta_boot <- lapply(fit_b$eva.boot, function(x) x$para[(1:n_x)])

    direct_est_boot2 <- list()
    for(i_boot in 1:nboot){
      .boot2 <- fit_b$eva_boot_direct[[i_boot]]
      .direct_est2 <- cbind(t(.boot2$metric), C = NA, sen = t(.boot2$senspe[,1]), spe = t(.boot2$senspe[,2]))
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
  est_em_like = c(b = fit_em_like$b, fit_em_like$beta, fit_em_like$eva),
  direct = c(b = simu$gamma, beta = simu$beta, direct_est),
  direct_b = c(b = simu$gamma, beta = simu$beta, direct_est_b)
  )

  if(var.smcure == T){
    b_sd0    <- apply(fit$b_boot, 2, function(x) c( mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975))) )
    beta_sd0 <- apply(fit$beta_boot, 2, function(x) c( mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975))) )

    boot1_res <- apply(direct_est_boot1, 2, function(x) c( boot1.mean = mean(x), boot1.sd = sd(x), boot1 = quantile(x, c(0.025, 0.975), na.rm = TRUE)) )
    boot2_res <- apply(direct_est_boot2, 2, function(x) c( boot2.mean = mean(x), boot2.sd = sd(x), boot2 = quantile(x, c(0.025, 0.975), na.rm = TRUE)) )

    res_sd <- rbind(
                 cbind(b_sd0, beta_sd0, fit$eva_sd),
                 boot1_res,
                 boot2_res)

  }

  # Evaluation at test dataset
  eva_test <- eva_cure(time = simu_test$surv[,1], delta = simu_test$surv[,2], X = simu_test$X,beta = fit$beta,
                       Z = simu_test$Z, b = fit$b, surv = fit$s, model = "PH", baseline = T)
  eva_test$sensep <- cutSenspe(eva_test$sensep, cutpoint )
  test_est <- c(eva_test$metric, sen = eva_test$sensep[,1], sep = eva_test$sensep[,2])

  eva_direct_test <- eva_cure_direct(time = simu_test$surv[,1], delta = simu_test$surv[,2], X = simu_test$X,beta = fit$beta,
                              Z = simu_test$Z, b = fit$b, surv = fit$s, model = "PH", baseline = T, cutpoint = cutpoint)
  test_direct <- cbind(t(eva_direct_test$metric), C = NA, sen = t(eva_direct_test$senspe[,1]), spe = t(eva_direct_test$senspe[,2]))

  eva_direct_test_b <- eva_cure_direct(time = simu_test$surv[,1], delta = simu_test$surv[,2], X = simu_test$X,beta = fit$beta,
                                     Z = simu_test$Z, b = fit$b, surv = fit$s, model = "PH", baseline = T, cutpoint = cutpoint, posterior = FALSE)
  test_direct_b <- cbind(t(eva_direct_test_b$metric), C = NA, sen = t(eva_direct_test_b$senspe[,1]), spe = t(eva_direct_test_b$senspe[,2]))


  res_test <-rbind(test_est = test_est, test_direct = test_direct, test_direct_b = test_direct_b)
  rownames(res_test) <- c("test_est", "test_direct", "test_direct_b")

  list(res = res, info = simu$info, res_test = res_test, res_sd = res_sd)
}

