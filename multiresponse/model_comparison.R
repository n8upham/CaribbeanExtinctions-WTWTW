library(brms)

##remove prior data
rm(list=ls()) 

##models with different response links
##which one is better response link?
##this yields warnings about using k_fold instead but leave as is
load("stan1multi_method.Rdata")
mod_cum <- add_criterion(mod_cum, c("loo", "waic"))
mod_aca <- add_criterion(mod_aca, c("loo", "waic"))
mod_cra <- add_criterion(mod_cra, c("loo", "waic"))
mod_sra <- add_criterion(mod_sra, c("loo", "waic"))
save.image("stan1multi_method.Rdata")

##print comparison 0 is best model
sink("loo_criterion_links.txt")
print(loo_compare(mod_cum, mod_aca, mod_cra, mod_sra, criterion = "loo"))
sink()

##remove prior data
rm(list=ls()) 

##models with different response variable
load("stan1multi_model.Rdata")
mod22 <- add_criterion(mod22, c("loo", "waic"))
mod24 <- add_criterion(mod24, c("loo", "waic"))
modf <- add_criterion(modf, c("loo", "waic"))
modi22 <- add_criterion(modi22, c("loo", "waic"))
modi24 <- add_criterion(modi24, c("loo", "waic"))
save.image("stan1multi_model.Rdata")

##print comparison 0 is best model
sink("loo_criterion_model_cumulative.txt")
print(loo_compare(mod22, mod24, criterion = "loo"))
print(loo_compare(modf,modi22,modi24, criterion = "loo"))
sink()

##remove prior data
rm(list=ls()) 

##models with different response variable
load("stan1multi_model_aca.Rdata")
mod22 <- add_criterion(mod22, c("loo", "waic"))
modi22 <- add_criterion(modi22, c("loo", "waic"))
save.image("stan1multi_model_aca.Rdata")

##print loo values, cannot be compared bc of different data sets
sink("loo_criterion_model_acat.txt")
print(mod22$criteria$loo)
print(modi22$criteria$loo)
sink()

##remove all data
rm(list=ls()) 