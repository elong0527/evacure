#' Semiparametric mixture cure model
#'
#' The function is same as smcure:::smcure
#' that fit a mxture cure models. This one output information for bootstrap
#' samples. See options in smcure for the detail.
#'
#' @export
smcure1 <-
  function(formula,cureform,offset=NULL,data,na.action=na.omit,model= c("aft", "ph"),link="logit", Var=TRUE,emmax=50,eps=1e-7,nboot=100,
           init = NULL, em = "smcure1", cutpoint = c(0.1,0.25,0.5,0.75,0.9), eva_model = NULL)
  {
    # if(em == "smcure") {
    #   print("smcure:::em")
    #   em <- smcure:::em }else{
    #     print("revised em")
    #   }

    if(model == "ph"){
      eva_model = "PH"
    }

    if(model == "aft" & is.null(eva_model) ){
      stop("eva_model is required to calculate prognostic accuracy for AFT modeld")
    }

    call <- match.call()
    model <- match.arg(model)
    cat("Program is running..be patient...")
    ## prepare data
    data <- na.action(data)
    n <- dim(data)[1]
    mf <- model.frame(formula,data)
    cvars <- all.vars(cureform)
    Z <- as.matrix(cbind(rep(1,n),data[,cvars]))
    colnames(Z) <- c("(Intercept)",cvars)
    if(!is.null(offset)) {
      offsetvar <- all.vars(offset)
      offsetvar<-data[,offsetvar]
    }else offsetvar <- NULL
    Y <- model.extract(mf,"response")
    X <- model.matrix(attr(mf,"terms"), mf)
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    Time <- Y[,1]
    Status <- Y[,2]
    bnm <- colnames(Z)
    nb <- ncol(Z)
    if(model == "ph") {
      betanm <- colnames(X)[-1]
      nbeta <- ncol(X)-1}
    if(model == "aft"){
      betanm <- colnames(X)
      nbeta <- ncol(X)}
    ## initial value
    w <- Status
    b <- eval(parse(text = paste("glm", "(", "w~Z[,-1]",",family = quasibinomial(link='", link, "'",")",")",sep = "")))$coef
    if(model=="ph") beta <- coxph(Surv(Time, Status)~X[,-1]+offset(log(w)), subset=w!=0, method="breslow")$coef
    if(model=="aft") beta <- survreg(Surv(Time,Status)~X[,-1])$coef
    ## do EM algo
    if(! is.null(init)){
      w <- init$w
      b <- init$b
      beta <- init$beta
    }
    emfit <- em(Time,Status,X,Z,offsetvar,b,beta,model,link,emmax,eps)
    b <- emfit$b
    beta <- emfit$latencyfit
    s <- emfit$Survival
    logistfit <- emfit$logistfit

    if(model == "ph"){
      eva <- eva_cure(Time, Status, X[, -1], beta, Z, b, s, model = "PH")
      eva$sensep <- cutSenspe(eva$sensep, cutpoint )
      eva <- c(eva$metric, sen = eva$sensep[,1], sep = eva$sensep[,2])
    }
    if(model == "aft"){
      eva <- eva_cure(Time, Status, X, - beta, Z, b, s, model = eva_model)
      eva$sensep <- cutSenspe(eva$sensep, cutpoint )
      eva <- c(eva$metric, sen = eva$sensep[,1], sep = eva$sensep[,2])
    }

    ## Bootstrap begin
    if(Var){
      if(model=="ph") {b_boot<-matrix(rep(0,nboot*nb), nrow=nboot)
                       beta_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)
                       iter <- matrix(rep(0,nboot),ncol=1)}

      if(model=="aft") {b_boot<-matrix(rep(0,nboot*nb), nrow=nboot)
                        beta_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)}
      tempdata <- cbind(Time,Status,X,Z)
      data1<-subset(tempdata,Status==1);data0<-subset(tempdata,Status==0)
      n1<-nrow(data1);n0<-nrow(data0)

      eva.boot <- list()
      i<-1
      while (i<=nboot){
        id1<-sample(1:n1,n1,replace=TRUE);id0<-sample(1:n0,n0,replace=TRUE)
        bootdata<-rbind(data1[id1,],data0[id0,])
        bootZ <- bootdata[,bnm]
        if(model=="ph") bootX <- as.matrix(cbind(rep(1,n),bootdata[,betanm]))
        if(model=="aft") bootX <- bootdata[,betanm]
        bootfit <- em(bootdata[,1],bootdata[,2],bootX,bootZ,offsetvar,b,beta,model,link,emmax,eps)
        b_boot[i,] <- bootfit$b
        beta_boot[i,] <- bootfit$latencyfit

        if(model=="ph"){
          eva.boot[[i]] <- eva_cure(bootdata[,1],bootdata[,2],bootX[,-1],bootfit$latencyfit,bootZ,bootfit$b, bootfit$Survival, model = "PH")
          eva.boot[[i]]$sensep <- cutSenspe(eva.boot[[i]]$sensep, cutpoint )
        }
        if(model=="aft"){
          eva.boot[[i]] <- eva_cure(bootdata[,1],bootdata[,2],bootX, - bootfit$latencyfit,bootZ,bootfit$b, bootfit$Survival, model = eva_model)
          eva.boot[[i]]$sensep <- cutSenspe(eva.boot[[i]]$sensep, cutpoint )
        }

        if (bootfit$tau<eps) i<-i+1}
      b_var <- apply(b_boot, 2, var)
      beta_var <- apply(beta_boot, 2, var)
      b_sd <- apply(b_boot, 2, function(x) c( mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975))) )
      beta_sd <- apply(beta_boot, 2, function(x) c( mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975))) )



      res_boot <- do.call(rbind, lapply(eva.boot, function(x)
        tmp <- c(x$metric, sen = x$sensep[,1], sep = x$sensep[,2]) ) )
      eva_sd <- apply(res_boot, 2, function(x) c( mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.975))) )

    }
    fit<-list()
    class(fit) <- c("smcure")
    fit$logistfit <- logistfit
    fit$b <- b
    fit$beta <- beta
    fit$eva    <- eva
    fit$model <- model
    fit$status <- Status
#     browser()
    if(model == "ph"){fit$X <- X[, -1]}else { fit$X <- X }
    fit$X <- as.matrix(fit$X)
    fit$Z <- Z

    if(Var){
      fit$b_var <- b_var
      fit$b_sd <- b_sd
      fit$b_zvalue <- fit$b/b_sd
      fit$b_pvalue <- (1-pnorm(abs(fit$b_zvalue)))*2
      fit$beta_var <- beta_var
      fit$beta_sd <- beta_sd
      fit$beta_zvalue <- fit$beta/beta_sd
      fit$beta_pvalue <- (1-pnorm(abs(fit$beta_zvalue)))*2
      fit$b_boot <- b_boot
      fit$beta_boot <- beta_boot
      fit$eva_sd <- eva_sd
      fit$eva.boot <- eva.boot
    }
    cat(" done.\n")
    fit$call <- call
    fit$bnm <- bnm
    fit$betanm <- betanm
    fit$s <- s
    fit$Time <- Time
    if(model=="aft"){
      error <- drop(log(Time)-beta%*%t(X))
      fit$error <- error}
    if(model == "ph"){
      eva <- eva_smcure(fit, model = "PH")
      fit$metric <- eva$metric
      fit$sensep <- eva$sensep
    }
    fit
#     printsmcure(fit,Var)
  }

#' Rank function
#'
#' Rank estimating equation used in the M-step of the EM algorithm for the AFT mixture cure model.
#'
#' @param beta unknown parameters corresponding to latency part
#' @param Time time to event of interest
#' @param X a vector or matrix of covariates corresponding to latency part
#' @param n total number of observations
#' @param w conditional probability of the individual remaining uncured
#' @param Status censoring indicator, 1=event of interest happens, and 0=censoring
smrank <- function (beta, Time, X, n, w, Status)
{
  error <- drop(log(Time) - beta %*% t(X))
  tp <- numeric()
  for (i in 1:n) {
    tp[i] <- sum(as.numeric((error[i] - error) < 0) * abs(error[i] -
                                                            error) * w * Status[i])
  }
  sum(tp)/n
}

#' Evaluate Mixture cure model fitted by smcure
#'
#' @export
#' @param time survival time
#' @param delta event status (1 = event)
#' @param X design matrix of survival part
#' @param beta parameters of survival part
#' @param Z design matrix of cure part
#' @param b parameters of cure part
#' @param model the type of survival model ("PH", "PO","Normal")
#' @param baseline whether the surv is baseline survival probability
eva_cure <- function(time,delta,X,beta,Z,b, surv, model, baseline = T){
  est.risk <- X %*% beta
  est.odds <- Z %*% b
  est.pi   <- logit.inv(est.odds)
  if(model == "PH" & baseline == T){ est.surv <- surv ^ exp( est.risk) }else{ est.surv <- surv}
  est.w <- w.cure(est.pi, delta, est.surv)

  k1  <- k.ind(risk = est.risk, pi = est.pi, model)
  k1w <- k.ind(risk = est.risk, pi = est.w, model)
  c1  <- cindex(risk = est.risk, pi = est.pi, time, delta )
  c1w <- cindex(risk = est.risk, pi = est.w,  time, delta )
  a1  <- auc.cure(est.odds, est.risk, est.surv, delta)
  res1 <- sen.spe.cure(est.odds, est.risk, est.surv, delta, verbose = F )

  res <- list( para = c( beta, b),
               metric = c(AUC = a1, K = k1, C = c1w),
               sensep = res1 )
  res
}

#' Print smcure object with evaluation information
#'
#' Refer printsmcure() for detail information
#' @export
printsmcure1 <- function(x, Var = TRUE, ROC = TRUE, ...)
{
  if (is.null(Var))
    Var = TRUE
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nCure probability model:\n")
  if (Var) {
    b <- array(x$b, c(length(x$b), 4))
    rownames(b) <- x$bnm
    colnames(b) <- c("Estimate", "Std.Error", "Z value",
                     "Pr(>|Z|)")
    b[, 2] <- x$b_sd
    b[, 3] <- x$b_zvalue
    b[, 4] <- x$b_pvalue
  }
  if (!Var) {
    b <- array(x$b, c(length(x$b), 1))
    rownames(b) <- x$bnm
    colnames(b) <- "Estimate"
  }
  print(b)
  cat("\n")
  cat("\nFailure time distribution model:\n")
  if (Var) {
    beta <- array(x$beta, c(length(x$beta), 4))
    rownames(beta) <- x$betanm
    colnames(beta) <- c("Estimate", "Std.Error", "Z value",
                        "Pr(>|Z|)")
    beta[, 2] <- x$beta_sd
    beta[, 3] <- x$beta_zvalue
    beta[, 4] <- x$beta_pvalue
  }
  if (!Var) {
    beta <- array(x$beta, c(length(x$beta), 1))
    rownames(beta) <- x$betanm
    colnames(beta) <- "Estimate"
  }
  print(beta)
  cat("\n")
  cat("\n Evaluation Metrics\n")
  if (!Var){
    metric <- x$metric
    if(ROC){
      plot(1-x$sensep[,2], x$sensep[,1], type = "s", xlab = "1 - Specificity", ylab = "Sensitivity")
    }
  }
  if (Var){
    ci <- apply(x$metric.boot, 2, quantile, probs = c(0.025, 0.975))
    metric <- cbind(x$metric, t(ci))
    colnames(metric)[1] <- c("Estimate")
    if(ROC){
      library(plyr)
      sensep.res <- lapply(x$sensep.boot,
                           function( sensep ){
                             cutpoints <- seq(0,1-1e-5,length = 100)
                             cutSenspe(sensep, cutpoints)
                           }
      )
      type <- rep(1:100, length(sensep.res))
      sensep.res <- do.call(rbind, sensep.res)
      rownames(sensep.res) <- NULL
      sensep.res <- data.frame(type, sensep.res)
      sen.ci <- ddply(sensep.res, .(type), function(res) c( quantile(res$sen, c(0.025, 0.975)), quantile(res$spe, 0.5 ) ) )

      plot(1-fit$sensep[,2], fit$sensep[,1], type = "s", xlab = "1 - Specificity", ylab = "Sensitivity")
      lines(1 - sen.ci[,4], sen.ci[,2], lty = 2)
      lines(1 - sen.ci[,4], sen.ci[,3], lty = 2)
    }
  }
  print(metric)
  invisible(x)
}

#' Evaluate Mixture cure model fitted by smcure
#'
#' @export
#' @param fit an extended smcure object
#' @param model the type of survival model ("PH", "PO","Normal")
eva_smcure <- function(fit, model){

  X <- fit$X
  Z <- fit$Z
  time <- fit$Time
  delta <- fit$status
  beta <- fit$beta
  b <- fit$b
  surv <- fit$s
  eva_cure(time,delta,X,beta,Z,b,surv,model)

}

