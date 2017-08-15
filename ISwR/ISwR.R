library(ISwR)     # Get melanoma dataset
library(survival)
library(evacure)

# Transfer covariates
data(melanom)
melanom$year <- melanom$days / 365.25
melanom$l.thick <- log(melanom$thick)
melanom$sex1    <- melanom$sex == 2 #( 1 = male)
melanom$ulc1     <- melanom$ulc == 1 #( 1 = present)

median(melanom$year) # median survival time
1 - sum(melanom$status == 1) / nrow(melanom) # Censoring proportion


##-----------------------------------------
## Kaplan-Meier Estimator Curve
##-----------------------------------------

# Ulceration
a <- survfit(Surv(melanom$days/365.25, melanom$status == 1) ~ ulc, data = melanom)
plot(a, ylab = "Survival Probability", xlab = "Year", lty = c(2,1), mark.time = T)
legend("bottomleft", lty = c(1,2), legend = c("No","Yes"))

# Tumor thickness
b <- survfit(Surv(melanom$days/365.25, melanom$status == 1) ~ I(l.thick > median(l.thick) ) , data = melanom)
plot(b, ylab = "Survival Probability", xlab = "Year", lty = c(2,1), mark.time = T)
legend("bottomleft", lty = c(1,2), legend = c("Above Median","Below Median"))

##-----------------------------------------
## CoxPH cure model with K-index
##-----------------------------------------
fit <- smcure1( Surv(melanom$days, melanom$status == 1) ~ sex1 + l.thick + ulc1,
               cureform = ~ sex1 + l.thick + ulc1, data = melanom, model = "ph",
               Var = T, em = "smcure"
)
printsmcure1(fit, Var = T)

# Fitted values
rbind( est = fit$beta, fit$beta_sd[-1, ])  # Cox part
rbind( est = fit$b, fit$b_sd[-1,] )        # Logistic part

# K index
c(fit$eva[2] , fit$eva_sd[-1,2])
