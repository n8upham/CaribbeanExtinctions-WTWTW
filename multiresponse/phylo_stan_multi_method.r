library(MASS)
library(geiger)
library(MCMCglmm)
library(MCMCvis)
library(brms)
library(parallel)
##using stan instead of phyr/inla for multinomial response
##see http://www.r-inla.org/faq#TOC-I-have-multinomial-data-but-I-cannot-find-the-multinomial-likelihood-isn-t-it-available-

##remove prior data
rm(list=ls()) 

##use the trees
load("ultra_trees.Rdata")

##read & clean up data
dat<-read.csv("Final_WTWTW_dataset_July2020.csv")

##reconfigure dat1
dat1<-dat[,]

##scale and transform data
##response and species vars
dat1$y<-dat1$Graded_extinction_response+1
dat1$mass<-scale(log10(dat1$Body_mass))
dat1$mas2<-dat1$mass^2

##physical vars
dat1$area<-scale(log10(dat1$Area_km2))
dat1$elev<-scale(log10(dat1$Mean_Elevation+1))
dat1$maxe<-scale(log10(dat1$Max_Elevation)) 
dat1$volc<-dat1$Active_Holocene_volcano ##this is binary 
dat1$hurr<-scale(dat1$Hurricane) 

##human vars
dat1$fore<-scale(dat1$Proportion_ForestCover2000)
dat1$loss<-scale(dat1$Proportion_ForestLoss2000)
dat1$impa<-scale(dat1$Mean_HFP)
dat1$fhar<-scale(dat1$First_human_arrival)
dat1$mong<-dat1$Mongoose ##this is binary

##must "as.character" to match the species name to tree labels in phyr
##otherwise the species names and trees don't match
dat1$sp = as.character(dat1$taxon)

##there is code for running across all trees elsewhere
##run  for only 1 tree
##this works like a charm compared to past scripts
A <- ape::vcv.phylo(phylo[[50]])

##new brms notation requires gr structure
##leave only some columns to reduce missingness
dat1<-dat1[,c('y', 'mass', 'mas2', 'area', 'elev', 'maxe', 'fore', 'loss', 'impa', 'fhar', 'mong', 'sp', 'Island')]

##somehow these got transformed to matrices
##make all numeric
dat1$mass<-as.numeric(dat1$mass)
dat1$mas2<-as.numeric(dat1$mas2)
dat1$area<-as.numeric(dat1$area)
dat1$elev<-as.numeric(dat1$elev)
dat1$maxe<-as.numeric(dat1$maxe)
dat1$fore<-as.numeric(dat1$fore)
dat1$loss<-as.numeric(dat1$loss)
dat1$impa<-as.numeric(dat1$impa)
dat1$fhar<-as.numeric(dat1$fhar)

##no imputation

##bayesian computation
##both species & island fixed effects
##use only one tree for now
##include 2/2 vars species/island
##cumulative
mod_cum<-brm(y ~ mass + mas2 + area + elev + (1|gr(sp, cov = A)) + (1|Island), data = dat1, data2 = list(A = A), family = "cumulative", cores=numCores, iter=10000, chains=2, save_all_pars=T, control = list(adapt_delta = 0.9))

sink("multi_stan_response.txt")
print("mod_cum")
print(summary(mod_cum))
sink()

##save the models
save.image("stan1multi_method.Rdata")

##adjacent categorical
mod_aca<-brm(y ~ mass + mas2 + area + elev + (1 | gr(sp, cov = A))  + (1|Island), data = dat1, data2 = list(A = A), family = "acat", cores=numCores, iter=10000, chains=2, save_all_pars=T, control = list(adapt_delta = 0.9))

sink("multi_stan_response.txt", append=T)
print("mod_aca")
print(summary(mod_aca))
sink()

##save the models
save.image("stan1multi_method.Rdata")

##continuation ratio
mod_cra<-brm(y ~ mass + mas2 + area + elev + (1 | gr(sp, cov = A))  + (1|Island), data = dat1, data2 = list(A = A), family = "cratio", cores=numCores, iter=10000, chains=2, save_all_pars=T, control = list(adapt_delta = 0.95))

sink("multi_stan_response.txt", append=T)
print("mod_cra")
print(summary(mod_cra))
sink()

##save the models
save.image("stan1multi_method.Rdata")

##stopping ratio
mod_sra<-brm(y ~ mass + mas2 + area + elev + (1 | gr(sp, cov = A))  + (1|Island), data = dat1, data2 = list(A = A), family = "sratio", cores=numCores, iter=20000, chains=2, save_all_pars=T, control = list(adapt_delta = 0.9))

sink("multi_stan_response.txt", append=T)
print("mod_sra")
print(summary(mod_sra))
sink()

##save the models
save.image("stan1multi_method.Rdata")

##remove all the data
rm(list=ls()) 