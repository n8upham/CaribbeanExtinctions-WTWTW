library(MASS)
library(geiger)
library(MCMCglmm)
library(MCMCvis)
library(phyr)
library(INLA)
library(parallel)
library(ppcor)
library(reshape2)

##model has random intercepts for island and species
##b series has an intercept for whole regression
##Thank Daijiang Li (phyr) for figuring out how to make model converge
##Thank Russell Dinnage for the inla implementation after phyr got rid of the random effects

##remove prior data
rm(list=ls())

##load data
load("phyr100binomial_dat.Rdata") 

##parallelize
##detect cores
numCores <- detectCores()

### mclapply regressions
##define function
##pglmm uses notation __ for phylogenetic effects, bayes=T uses inla
##significant full model
f22 <- function(i) {
  pglmm(sigf, data = dat1, cov_ranef=list(sp = phylo[[i]]), family="binomial", bayes = T, marginal.summ="median", calc.DIC=T, prior= "pc.prior.auto")
}

##run and time the function 
system.time({
##must load the library to run on the forks
  library(phyr)
  save22 <- mclapply(1:100, f22, mc.cores =numCores)
})

##save the models
save.image("phyr100binomial_v3.Rdata")

##another regression
##two more vars for island
##significant imputed model 
f24 <- function(i) {
  pglmm(sigi, data = dat2, cov_ranef=list(sp = phylo[[i]]), family="binomial", bayes = T, marginal.summ="median", calc.DIC=T, prior= "pc.prior.auto")
}

system.time({
  library(phyr)
  save24 <- mclapply(1:100, f24, mc.cores =numCores)
})

##save the models
save.image("phyr100binomial_v3.Rdata")

##another regression
##all variables
fma <- function(i) {
  pglmm(full, data = dat1, cov_ranef=list(sp = phylo[[i]]), family="binomial", bayes = T, marginal.summ="median", calc.DIC=T, prior= "pc.prior.auto")
}

system.time({
  library(phyr)
  savema <- mclapply(1:100, fma, mc.cores =numCores)
})

##save the models
save.image("phyr100binomial_v3.Rdata")
