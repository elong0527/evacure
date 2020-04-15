#' Evaluate the direct estimator
#'
#' @export
eva_cure_direct <- function(time, delta, X, beta, Z, b, surv, model, cutpoint, n_post = 500, baseline = TRUE){
  n <- length(time)

  est.risk <- X %*% beta
  est.odds <- Z %*% b
  est.pi   <- logit.inv(est.odds)
  if(toupper(model) == "PH" & baseline == T){ est.surv <- surv ^ exp( est.risk) }else{ est.surv <- surv}
  est.w <- w.cure(est.pi, delta, est.surv)

  .eva_cure_observed <- function(risk, odds, uncure, model, cutpoint){

    k <- k.ind(risk = risk[uncure == 1, ], pi = rep(1, sum(uncure)), model = model)

    auc <- performance( prediction(odds, uncure ), measure = "auc")@y.values[[1]]
    sen <- performance( prediction(odds, uncure ),  "sens" )
    spe <- performance( prediction(odds, uncure ),  "spec" )
    res <- cbind(sen = sen@y.values[[1]], spe = spe@y.values[[1]], cut = spe@x.values[[1]])
    res1<- cutSenspe(res, cutpoint )
    res1 <- data.frame(res1)
    res1$cut <- rownames(res1)
    list(k = k, auc = auc,
         senspe = res1)
  }

  res_post <- list()
  for(i_post in 1:n_post){

    y_imp    <- ifelse(delta == 1, 1, rbinom(n, size = 1, prob = est.w)) # Uncure imputation

    .res <- .eva_cure_observed(est.risk, est.odds, uncure = y_imp, model = model, cutpoint = cutpoint)
    res_post[[i_post]] <- .res
  }


  post0_metric <- do.call(rbind, lapply(res_post, function(x) c(AUC = x$auc, k = x$k)))
  post0_senspe <- do.call(rbind, lapply(res_post, function(x) x$senspe))

  post_metric <- apply(post0_metric, 2, mean)
  post_senspe <- data.frame(sen = tapply(post0_senspe[,1], post0_senspe$cut, function(x) mean(x)),
                            spe = tapply(post0_senspe[,2], post0_senspe$cut, function(x) mean(x)))
  post_senspe$cut <- rownames(post_senspe)

  post_fit <- list(metric = post_metric,
                   senspe = post_senspe)

  post_fit

}
