#---------------------
# Simulation for parametric Cox PH model
#---------------------

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

# Set up Simulation Environment

# task_id <- 1
set.seed(task_id * 5)

#-----------------------------------
library(devtools)
# library(evacure)
#-----------------------------------

devtools::load_all()
source("simu_data.R")       ## Simulate data
source("simu.coxph.R") ## Run simulation with Cox PH model

# set.seed(1234)
niter <- 1            ## Number of simulation replication
var.smcure = TRUE
nboot = 200
n_post = 500

res1 = list()
res2 = res3 = res1
for( iter in 1:niter){

  print(iter)

  ## 60/50
  res1[[iter]] <- simu.coxph(N = 200, c.min = 0, c.max = 18, model ="PH",
                     .beta = c(-1, 1), .gamma = c(0.5, 1, -1), share = T,
                     var.smcure = var.smcure, nboot = nboot, n_post = n_post
  )

  ## 60/30
  res2[[iter]] <- simu.coxph(N = 200, c.min = 0, c.max = 3, model ="PH",
                     .beta = c(-1, 1), .gamma = c(2, 1, -1), share = T,
                     var.smcure = var.smcure, nboot = nboot, n_post = n_post
  )

  ## 80/50
  res3[[iter]] <- simu.coxph(N = 200, c.min = 0, c.max = 2, model ="PH",
                     .beta = c(-1, 1), .gamma = c(0.5, 1, -1), share = T,
                     var.smcure = var.smcure, nboot = nboot, n_post = n_post
  )


}

filename <- paste0(task_id,".Rdata")
print(filename)
save(res1, res2, res3, file = filename)

#
# library(purrr)
# library(dplyr)
# library(tidyr)
#
# path <- "/SFS/scratch/zhanyilo/LDA"
#
#
# res_1 <- list()
# for(i in 1:1000){
#   try({
#     load(file.path(path, paste0(i, ".Rdata")))
#     res_1[[i]] <- res1[[1]]
#   })
# }
#
#
# db <- res_1
# tmp <- do.call(rbind, map(db, "res"))
# type <- rownames(tmp)
# tmp <- data.frame(type = type, tmp)
#
# # MSE
# tmp_long <- tmp %>% subset(type %in% c("true", "est", "direct", "boot1.mean", "boot2.mean")) %>%
#   pivot_longer(cols = b1:`spe.90.`)
#
#
# true_value <- tmp_long %>% subset(type == "true") %>%
#   group_by(name) %>%
#   summarise(true = mean(value))
#
# est_sd <-   tmp %>% group_by(type) %>% summarise_if(is.numeric, mean) %>%
#   subset(type %in% c("sd", "boot1.sd", "boot2.sd")) %>%
#   select(type, AUC, sen.25., sen.50., sen.75., spe.25., spe.50., spe.75.) %>%
#   mutate(type = factor(type, levels = c("sd", "boot1.sd", "boot2.sd"),
#                        labels = c("est", "boot1.mean", "boot2.mean"))) %>%
#   pivot_longer(cols = AUC:spe.75., values_to = "sd_est") %>%
#   arrange(type)
#
#
# tmp_res <- tmp_long %>% left_join(true_value) %>%
#   group_by(type, name) %>%
#   summarise(true = mean(true),
#             mean = mean(value),
#             sd_empirical = sd(value),
#             rmse = sqrt(mean( (value - true)^2 ))) %>%
#   subset(name %in% c("AUC", "sen.25.", "sen.50.", "sen.75.", "spe.25.", "spe.50.", "spe.75.")) %>%
#   ungroup() %>%
#   left_join(est_sd) %>%
#   mutate(type = factor(type, levels = c("true", "est", "direct", "boot1.mean", "boot2.mean"))) %>%
#   subset(type != "true") %>%
#   arrange(name, type)
#
# tmp_res
