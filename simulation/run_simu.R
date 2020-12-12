#---------------------
# Simulation for parametric Cox PH model
#---------------------

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

# Set up Simulation Environment

# task_id <- 1
set.seed(task_id * 5)

#-----------------------------------
library(devtools)
library(evacure)
#-----------------------------------

# devtools::load_all()
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

#----------------------
# HPC Submission code
#----------------------

# cd /SFS/scratch/zhanyilo/boundary
# rm *
# cp ~/evacure/simulation/*.R .
# module add R/3.6.1
# qsub -t 1:2000 ~/runr.sh run_simu.R


#-------------------------------
# Summarize simulation results
#-------------------------------
# library(purrr)
# library(dplyr)
# library(tidyr)
#
# path <- "/SFS/scratch/zhanyilo/evacure"
#
#
# res_1 <- list()
# res_2 <- list()
# res_3 <- list()
# for(i in 1:1000){
#   try({
#     load(file.path(path, paste0(i, ".Rdata")))
#     res_1[[i]] <- res1[[1]]
#     res_2[[i]] <- res2[[1]]
#     res_3[[i]] <- res3[[1]]
#   })
# }
#
#
# db <- res_2
# tmp <- do.call(rbind, map(db, "res"))
# type <- rownames(tmp)
# tmp <- data.frame(type = type, tmp)
#
# tmp_sd <- do.call(rbind, map(db, "res_sd"))
# type <- rownames(tmp_sd)
# tmp_sd <- data.frame(type = type, tmp_sd)
# names(tmp_sd) <- gsub("sep", "spe", names(tmp_sd))
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
# # TEST data
# tmp_test <- do.call(rbind, map(db, "res_test"))
# type <- rownames(tmp_test)
# tmp_test <- data.frame(type = type, tmp_test)
# names(tmp_test) <- gsub("sep", "spe", names(tmp_test))
#
# est_test <-   tmp_test %>% group_by(type) %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>%
#   subset(type %in% c("test_est", "test_direct")) %>%
#   select(type, AUC, sen.25., sen.50., sen.75., spe.25., spe.50., spe.75.) %>%
#   mutate(type = factor(type, levels = c("test_est", "test_direct"),
#                        labels = c("est", "direct"))) %>%
#   pivot_longer(cols = AUC:spe.75., values_to = "est_test") %>%
#   arrange(type)
#
# # SD
# est_sd <-   tmp_sd %>% group_by(type) %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>%
#   subset(type %in% c("sd", "boot1.sd", "boot2.sd")) %>%
#   select(type, AUC, sen.25., sen.50., sen.75., spe.25., spe.50., spe.75.) %>%
#   mutate(type = factor(type, levels = c("sd", "boot1.sd", "boot2.sd"),
#                        labels = c("est", "direct", "boot2"))) %>%
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
#   left_join(est_test) %>%
#   mutate(type = factor(type, levels = c("true", "est", "direct", "boot1.mean", "boot2.mean"))) %>%
#   subset(type != "true") %>%
#   arrange(name, type)
#
# tmp_res %>% mutate_if(is.numeric, formatC, digits = 3, format = "f") %>% write.csv(file = "res2.csv")
#
# # Get simulation results
# t1 <- read.csv("res1.csv")
# t1$X <- 1
# t2 <- read.csv("res2.csv")
# t2$X <- 2
# t3 <- read.csv("res3.csv")
# t3$X <- 3
#
# bind_rows(t1, t2, t3) %>% select(X, type, name, true, mean, sd_empirical, sd_est, rmse, est_test) %>%
#   arrange(rev(type), X, name) %>% write.csv("simu_res.csv")
