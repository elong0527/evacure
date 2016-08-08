smsurv <- function (Time, Status, X, beta, w, model)
{
  death_point <- sort(unique(subset(Time, Status == 1)))
  if (model == "ph")
    coxexp <- exp((beta) %*% t(X[, -1]))
  lambda <- numeric()
  event <- numeric()
  for (i in 1:length(death_point)) {
    event[i] <- sum(Status * as.numeric(Time == death_point[i]))
    if (model == "ph")
      temp <- sum(as.numeric(Time >= death_point[i]) *
                    w * drop(coxexp))
    if (model == "aft")
      temp <- sum(as.numeric(Time >= death_point[i]) *
                    w)
    temp1 <- event[i]
    lambda[i] <- temp1/temp
  }
  HHazard <- numeric()
  for (i in 1:length(Time)) {
    HHazard[i] <- sum(as.numeric(Time[i] >= death_point) *
                        lambda)
    if (Time[i] > max(death_point))
#       HHazard[i] <- Inf
#       print("revised tail")
      HHazard[i] <- Time[i] #Exponential Tail
    if (Time[i] < min(death_point))
      HHazard[i] <- 0
  }
  survival <- exp(-HHazard)
  list(survival = survival)
}

em <- function (Time, Status, X, Z, offsetvar, b, beta, model, link,
                emmax, eps)
{
  w <- Status
  n <- length(Status)
  if (model == "ph")
    s <- smsurv(Time, Status, X, beta, w, model)$survival
  if (model == "aft") {
    if (!is.null(offsetvar))
      Time <- Time/exp(offsetvar)
    error <- drop(log(Time) - beta %*% t(X))
    s <- smsurv(error, Status, X, beta, w, model)$survival
  }
  convergence <- 1000
  i <- 1
  while (convergence > eps & i < emmax) {
    uncureprob <- matrix(exp((b) %*% t(Z))/(1 + exp((b) %*%
                                                      t(Z))), ncol = 1)
    if (model == "ph") {
      survival <- drop(s^(exp((beta) %*% t(X[, -1]))))
    }
    if (model == "aft") {
      error <- drop(log(Time) - beta %*% t(X))
      survival <- s
    }
    w <- Status + (1 - Status) * (uncureprob * survival)/((1 -
                                                             uncureprob) + uncureprob * survival)
    logistfit <- eval(parse(text = paste("glm", "(", "w~Z[,-1]",
                                         ",family = quasibinomial(link='", link, "'", ")",
                                         ")", sep = "")))
    update_cureb <- logistfit$coef
    if (!is.null(offsetvar))
      update_cureb <- as.numeric(eval(parse(text = paste("glm",
                                                         "(", "w~Z[,-1]+offset(offsetvar)", ",family = quasibinomial(link='",
                                                         link, "'", ")", ")", sep = "")))$coef)
    if (model == "ph") {
      update_beta <- coxph(Surv(Time, Status) ~ X[, -1] +
                             offset(log(w)), subset = w != 0, method = "breslow")$coef
      if (!is.null(offsetvar))
        update_beta <- coxph(Surv(Time, Status) ~ X[,
                                                    -1] + offset(offsetvar + log(w)), subset = w !=
                               0, method = "breslow")$coef
      update_s <- smsurv(Time, Status, X, beta, w, model)$survival
    }
    if (model == "aft") {
      update_beta <- optim(rep(0, ncol(X)), smrank, Time = Time,
                           X = X, n = n, w = w, Status = Status, method = "Nelder-Mead",
                           control = list(reltol = 1e-04, maxit = 500))$par
      update_s <- smsurv(error, Status, X, beta, w, model)$survival
    }
    convergence <- sum(c(update_cureb - b, update_beta -
                           beta)^2) + sum((s - update_s)^2)
    b <- update_cureb
    beta <- update_beta
    s <- update_s
    uncureprob <- matrix(exp((b) %*% t(Z))/(1 + exp((b) %*%
                                                      t(Z))), ncol = 1)
    i <- i + 1
  }
  em <- list(logistfit = logistfit, b = b, latencyfit = beta,
             Survival = s, Uncureprob = uncureprob, tau = convergence)
}
