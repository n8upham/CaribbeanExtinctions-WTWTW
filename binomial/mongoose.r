library(geiger)
library(phyr)
library(INLA)
library(reshape2)

##separated pre processing including first models from modeling across many trees
##this is the processing step

##model has random intercepts for island and species
##b series has an intercept for whole regression
##Thank Daijiang Li (phyr) for figuring out how to make model converge
##Thank Russell Dinnage for the inla implementation after phyr got rid of the random effects

##remove prior data
rm(list=ls()) 

##read & clean up data
dat<-read.csv("Final_WTWTW_dataset_July2020.csv")
dat1<-dat[,]

##scale and transform data
##response and species vars
dat1$y<-dat1$Binary_extinction_response
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
dat1$loss<-log10(dat1$Proportion_ForestLoss2000+1)
dat1$impa<-dat1$Mean_HFP
dat1$fhar<-dat1$First_human_arrival/1000
dat1$mong<-dat1$Mongoose ##this is binary


##must "as.character" to match the species name to tree labels in phyr
##otherwise the species names and trees don't match
dat1$sp = as.character(dat1$taxon)

##get phylogenies
phylo<-read.nexus("CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_resolvedUltra_July2020.nex")

##bayesian computation
##both species & island fixed effects

##full model, no imputation
##make formulae
##based on these results get rid of problematic variables maxe & mong
##switch mongoose for eleve
full<-as.formula(y ~ mass + mas2 + area + mong + volc + hurr + fhar + fore + loss + impa  + mass:fhar + (1|sp__) + (1|Island))

##use only one tree
modf<-pglmm(full, data = dat1, cov_ranef=list(sp = phylo[[50]]), family="binomial", bayes = T, marginal.summ="median", calc.DIC=T, prior= "pc.prior.auto")

##make a new dataset with imputed values
dat2<-dat1

## find missing values
mass.na <- which(is.na(dat1$mass))
hurr.na <- which(is.na(dat1$hurr))
impa.na <- which(is.na(dat1$impa))
fhar.na <- which(is.na(dat1$fhar))

##plug in imputation values from full model
dat2$mass[mass.na] <- modf$inla.model$summary.fitted.values[mass.na, "mean"]
dat2$mas2[mass.na] <- dat2$mass[mass.na]^2
dat2$hurr[hurr.na] <- modf$inla.model$summary.fitted.values[hurr.na, "mean"]
dat2$impa[impa.na] <- modf$inla.model$summary.fitted.values[impa.na, "mean"]
dat2$fhar[fhar.na] <- modf$inla.model$summary.fitted.values[fhar.na, "mean"]

##use only one tree imputed
modi<-pglmm(full, data = dat2, cov_ranef=list(sp = phylo[[50]]), family="binomial", bayes = T, marginal.summ="median", calc.DIC=T, prior= "pc.prior.auto")

##significant full model
sigf<-as.formula(y ~ mass + mas2 + elev + fhar + mass:fhar + (1|sp__) + (1|Island))

##significant imputed model
sigi<-as.formula(y ~ mass + mas2 + elev + fhar + mass:fhar + (1|sp__) + (1|Island))

##print out values
sink("mongoose_results_1_tree.txt")
print("no imputation")
print(summary(modf))
print("imputation")
print(summary(modi))
sink()

##save the models
save.image("phyr_mong_binomial_dat.Rdata")

