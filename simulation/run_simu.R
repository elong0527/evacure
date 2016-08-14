#---------------------
# Simulation for parametric Cox PH model
#---------------------

#-----------------------------------
library(devtools)
library(foreach)
library(plyr)

#-----------------------------------

library(evacure) 
source("simu_data.R")       ## Simulate data
source("simu.coxph.R") ## Run simulation with Cox PH model

set.seed(1234)         
niter <- 1             ## Number of simulation replication

res1 =list()
for( iter in 1:niter){
  ## Scenario 1
  ## censoring porpotion 60/  cure proportion 50
  res1[[iter]] <- simu.coxph(N = 200, c.min = 0, c.max = 1.2, model ="PH",
                     .beta = c(-1, 1), .gamma = c(0.5, 1, -1), share = T,
                     var.smcure = T
  )
  print(res1[[iter]]$res[, c(1:5,7)])  # Estiamted results 
  
}