#' Posterior of uncure probabilty
#'
#' @export
#' @param pi uncure probability
#' @param delta censoring status
#' @param surv survival probability
w.cure <- function(pi, delta, surv){
  delta + (1-delta) * logit.inv( log(surv) + logit(pi) )
#   delta + (1-delta) * pi * surv / (1 - pi + pi * surv)
}

#' Sensitivity and Specificity for mixture cure models
#'
#' @export
#' @param odds odds for cure part
#' @param risk risk for cox part
#' @param surv survival probability
#' @param delta censoring status
#' @param cut.point cut point for sensitivity (if NULL report all)
#' @param verbose  if TRUE, report AUC and draw ROC curve
sen.spe.cure <- function(odds, risk, surv, delta, cut.point = NULL, verbose = T){
  library(MASS)
  odds <- as.numeric(odds)
  risk <- as.numeric(risk)
  surv <- as.numeric(surv)
  delta <- as.numeric(delta)
  n <- length(odds)

  pi <- logit.inv(odds)
  w.est <- w.cure(pi, delta, surv)
  den <- sum(w.est)


  odds.o  <- sort(odds)
  w.est.o <- w.est[order(odds)]

  sen <- c(1, 1 - cumsum(w.est.o) / sum(w.est), 0)
  spe <- c(0,  cumsum(1 - w.est.o) / sum(1-w.est), 1 )
  odds.o <- c( min(odds.o) - 1e-5, odds.o, max(odds.o) + 1e-5)

  if(max(table(odds.o)) > 1){  # Handle ties
    if(verbose){ cat("there are ties in estimated odds\n") }
    sen <- unlist( tapply(sen, odds.o, function(x) x[c(length(x))]) )
    spe <- unlist( tapply(spe, odds.o, function(x) x[c(length(x))]) )
    a0 <- sum(diff(spe) * sen[-1] )
    a1 <- sum(diff(spe) * sen[-length(sen)] )
    auc0 <- mean( c(a0,a1) )
    sen[length(sen) - 1] <- sen[length(sen) - 2]
    sen[1] <- sen[2]
    sen <- c(1, sen)
    spe <- c(0, spe)

  }else{
    auc0 <- sum(diff(spe) * sen[-n] )
  }
  auc1 <- auc.cure(odds, risk, surv, delta) # Check the results is equivalent to directly caculation

  den <- ( sum(w.est) * sum(1-w.est) )
  rate <- ( den - sum(w.est * (1-w.est)) ) / den
  auc2 <- auc0 * rate

  if(verbose){
    cat(rate, "AUC is", round(auc0,3), "for integration with correction", round(auc2,3), ". AUC is ", round(auc1,3), "for directly caculation\n" )
    roc <- stepfun(rev( 1 - spe), c(rev(sen),1) )
    plot(roc, xlim = c(0,1), xlab = "1 - Specificity", ylab = "Sensitivity", main = "ROC curve" )
  }

  options(warn=-1)
  res <- cbind(sen, spe, cut = c(min(odds.o) - 1, unique(odds.o) ) )
  options(warn=0)

  if(! is.null(cut.point) ){
    res <- cut.senspe(res, cut.point)
  }
  res
}

#' Sensitivity and specificity for a cut point
#'
#' @export
#' @param results of sen.spe.cure() when cut.point = NULL
#' @param cut.point cutpoint for sensitivity and specificity
cutSenspe <- function(res, cut.point){
  cut.point <- quantile(res[,"cut"], cut.point)
  res <- res[order(res[,"cut"]),]
  res <- sapply(cut.point, function(cut.point){
    ind <- which(res[,"cut"] > cut.point)[1]
    colMeans( res[ind+ c(-1,0), ] )
  })
  t(res)
}

#' AUC for mixture cure models
#'
#' @param odds odds for cure part
#' @param risk risk for cox part
#' @param surv survival probability
#' @param delta censoring status
#' @export
auc.cure <- function(odds, risk, surv, delta){
  odds <- as.numeric(odds)
  risk <- as.numeric(risk)
  surv <- as.numeric(surv)
  delta <- as.numeric(delta)

  pi <- logit.inv(odds)
  w.est <- w.cure(pi, delta, surv)
  r1 <- outer(odds, odds, ">") + 0.5 * outer(odds, odds, "==")
  w1 <- outer(w.est, w.est, function(x,y) x * (1-y))
  num <- r1 * w1
  diag(num) <- 0
  den <- w1
  diag(den) <- 0

  sum(num) / sum(den)
}


