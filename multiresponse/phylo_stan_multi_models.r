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

##generating models with different predictor sets
##remove prior data
rm(list=ls()) 

##get ultrametric treess
load("ultra_trees.Rdata")

##get data already processed across methods
load("stan1multi_method.Rdata")

##clear heavy models from memory
rm("mod_cum", "mod_aca", "mod_cra", "mod_sra")

##no imputation

##bayesian computation
##both species & island fixed effects
##use only one tree for now
##include 2/2 vars species/island
mod22<-brm(y ~ mass + mas2 + area + elev + (1 | gr(sp, cov = A)) + (1|Island), data = dat1, data2 = list(A = A), family = "cumulative", cores=numCores, iter=10000)

sink("multi_stan_models.txt")
print("mod22")
print(summary(mod22))
sink()

##save the models
save.image("stan1multi_model.Rdata")

##include 2/4 vars species/island
mod24<-brm(y ~ mass + mas2 + area + elev + fore + mong + (1 | gr(sp, cov = A)) + (1|Island), data = dat1, data2 = list(A = A),  family = "cumulative", cores=numCores, iter=20000, control = list(adapt_delta = 0.95))

sink("multi_stan_models.txt", append=T)
print("mod24")
print(summary(mod24))
sink()

##save the models
save.image("stan1multi_model.Rdata")

##impute data 
##caution to the wind?
imp <- mice(dat1, m = 5, print = FALSE)

##for data2 must create a list of named lists of lenth equal to sets imputed from mice
Ai<-list()
for(i in 1:5){Ai[[i]]<-list(A=A)}

##use imputed data
##take out "+ impa + fhar" as their behavior suggests imputation has got rid of signal
##maxe gets chucked bc of collinearity

##parallelize across data sets
plan(multiprocess)

modf<-brm_multiple(y ~ mass + mas2 + area + elev + fore + loss + mong + (1 | gr(sp, cov = A)) + (1|Island), data = imp, data2 = Ai, family = "cumulative", cores=numCores, iter=50000, chains=2, control = list(adapt_delta = 0.95))

sink("multi_stan_models.txt", append=T)
print("modf")
print(summary(modf))
sink()

##save the models
save.image("stan1multi_model.Rdata")

##take out "+ impa + fhar" as their behavior suggests imputation has got rid of signal
##2/2 vars, imputed

##parallelize across data sets
plan(multiprocess)

modi22<-brm_multiple(y ~ mass + mas2 + area + elev + (1 | gr(sp, cov = A)) + (1|Island), data = imp, data2 = Ai, family = "cumulative", cores=numCores, iter=15000, chains=2, control = list(adapt_delta = 0.95))

sink("multi_stan_models.txt", append=T)
print("modi22")
print(summary(modi22))
sink()

##save the models
save.image("stan1multi_model.Rdata")

##take out "+ impa + fhar" as their behavior suggests imputation has got rid of signal
##2/4 vars, imputed

##parallelize across data sets
plan(multiprocess)

modi24<-brm_multiple(y ~ mass + mas2 + area + elev + fore + mong + (1 | gr(sp, cov = A)) + (1|Island), data = imp, data2 = Ai, family = "cumulative",cores=numCores, iter=25000, chains=2, control = list(adapt_delta = 0.95))

sink("multi_stan_models.txt", append=T)
print("modi24")
print(summary(modi24))
sink()

##save the models
save.image("stan1multi_model.Rdata")

##remove prior data
rm(list=ls()) 