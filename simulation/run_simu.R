#---------------------
# Simulation for parametric Cox PH model
#---------------------

i.par <- commandArgs(T)
i.par <- as.numeric(i.par)
root.file <- paste0(i.par)

#-----------------------------------
library(devtools)

#-----------------------------------

devtools::load_all()
source("simu_data.R")       ## Simulate data
source("simu.coxph.R") ## Run simulation with Cox PH model

# set.seed(1234)
niter <- 1            ## Number of simulation replication
var.smcure = TRUE
nboot = 200

res1 = list()
res2 = res3 = res1
for( iter in 1:niter){

  print(iter)

  ## 60/50
  res1[[iter]] <- simu.coxph(N = 200, c.min = 0, c.max = 18, model ="PH",
                     .beta = c(-1, 1), .gamma = c(0.5, 1, -1), share = T,
                     var.smcure = var.smcure, nboot = nboot
  )

  ## 60/30
  res2[[iter]] <- simu.coxph(N = 200, c.min = 0, c.max = 3, model ="PH",
                     .beta = c(-1, 1), .gamma = c(2, 1, -1), share = T,
                     var.smcure = var.smcure, nboot = nboot
  )

  ## 80/50
  res3[[iter]] <- simu.coxph(N = 200, c.min = 0, c.max = 2, model ="PH",
                     .beta = c(-1, 1), .gamma = c(0.5, 1, -1), share = T,
                     var.smcure = var.smcure, nboot = nboot
  )


}

save(res1, res2, res3, file = paste0(root.file,".Rdata") )

#
# library(purrr)
# library(dplyr)
# library(tidyr)
#
# summary_direct <- function(db){
#   t1 <- do.call(rbind, map(db, "res") )
#   t2 <- do.call(rbind, map(db, "info"))
#
#   type <- rownames(t1)
#   t1 <- data.frame(t1)
#   t1$type <- type
#
#   t1_long <- pivot_longer(t1,  cols = b1:`spe.90.`)
#
#   true_value <- t1_long %>% subset(type == "true") %>%
#     group_by(name) %>%
#     summarise(true = mean(value))
#
#   t1_res <- t1_long %>% left_join(true_value) %>%
#     subset(type %in% c("est", "direct")) %>%
#     group_by(type, name) %>%
#     summarise(true = mean(true),
#               mean = mean(value),
#               sd = sd(value),
#               rmse = sqrt(mean( (value - true)^2 ))) %>%
#     subset(name %in% c("AUC", "sen.25.", "sen.50.", "sen.75.", "spe.25.", "spe.50.", "spe.75.")) %>%
#     arrange(name, type)
#
#   t2_res <- apply(t2, 2, mean)
#
#   list(t1 = t1_res, t2 = t2_res)
# }
#
# round(res[c("true","est","direct","mean", "boot1.mean", "boot2.mean"),], 2)
# round(res[c("sd","boot1.sd","boot2.sd"),], 2)
#
# summary_direct(res1)
# summary_direct(res2)
# summary_direct(res3)
#
# res_x <- res1 %>% summary_direct()
# res_x$t1 %>% mutate_if(is.numeric, function(x) formatC(x, format = "f", digits = 3)) %>%
#   write.csv("res1.csv")
#
# res_x <- res2 %>% summary_direct()
# res_x$t1 %>% mutate_if(is.numeric, function(x) formatC(x, format = "f", digits = 3)) %>%
#   write.csv("res2.csv")
#
# res_x <- res3 %>% summary_direct()
# res_x$t1 %>% mutate_if(is.numeric, function(x) formatC(x, format = "f", digits = 3)) %>%
#   write.csv("res3.csv")
#
