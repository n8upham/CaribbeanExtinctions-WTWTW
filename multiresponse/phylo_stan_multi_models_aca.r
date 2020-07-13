library(MASS)
library(geiger)
library(MCMCglmm)
library(MCMCvis)
library(brms)
library(parallel)
library(mice)
library(future)
##using stan instead of phyr/inla for multinomial response
##see http://www.r-inla.org/faq#TOC-I-have-multinomial-data-but-I-cannot-find-the-multinomial-likelihood-isn-t-it-available-

##aca only runs with very restricted data

##remove prior data
rm(list=ls()) 

##get data already processed across methods
load("stan1multi_method.Rdata")

##clear heavy models from memory
rm("mod_cum", "mod_aca", "mod_cra", "mod_sra")

##no imputation

##bayesian computation
##both species & island fixed effects
##use only one tree for now
##include 2/2 vars species/island
mod22<-brm(y ~ mass + mas2 + area + elev + (1 | gr(sp, cov = A)) + (1|Island), data = dat1, data2 = list(A = A), family = "acat", cores=numCores, iter=10000, control = list(adapt_delta = 0.95))

sink("multi_stan_models_aca.txt")
print("mod22")
print(summary(mod22))
sink()

##save the models
save.image("stan1multi_model_aca.Rdata")

##impute data 
##caution to the wind?
imp <- mice(dat1, m = 5, print = FALSE)

##for data2 must create a list of named lists of lenth equal to sets imputed from mice
Ai<-list()
for(i in 1:5){Ai[[i]]<-list(A=A)}

##use imputed data
##take out "+ impa + fhar" as their behavior suggests imputation has got rid of signal
##maxe gets chucked bc of collinearity

##take out "+ impa + fhar" as their behavior suggests imputation has got rid of signal
##2/2 vars, imputed

##parallelize across data sets
plan(multiprocess)

modi22<-brm_multiple(y ~ mass + mas2 + area + elev + (1 | gr(sp, cov = A)) + (1|Island), data = imp, data2 = Ai, family = "acat", cores=numCores, iter=50000, chains=2, control = list(adapt_delta = 0.99))

sink("multi_stan_models_aca.txt", append=T)
print("modi22")
print(summary(modi22))
sink()

##save the models
save.image("stan1multi_model_aca.Rdata")

##remove prior data
rm(list=ls()) 