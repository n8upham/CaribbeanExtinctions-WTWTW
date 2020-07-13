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
dat1$loss<-scale(dat1$Proportion_ForestLoss2000)
dat1$impa<-scale(dat1$Mean_HFP)
dat1$fhar<-scale(dat1$First_human_arrival)
dat1$mong<-dat1$Mongoose ##this is binary

##test for multicollinearity
##missing value rows excluded

## focus only on island data characteristics
dat3<-dat1[!duplicated(dat1[,'Island']),]
mcol<- pcor(na.omit(dat3[,27:36]), method = "pearson")

##get corr coeff
m1<-mcol$estimate
diag(m1)<-NA	
m1.m<-melt(m1)

##get significantly correlated vars
m2<-mcol$p.value
diag(m2)<-NA	
m2.m<-melt(m2, value.name = "p.value")
m2.m<-subset(m2.m, p.value<0.05)

##add corr coeff
m2.m$corr<-m1.m[rownames(m2.m),'value']

##order by corr
m2.m<-m2.m[order(m2.m$corr),]

##export for decision making
write.csv(m2.m, "collinear.csv", row.names =F)

##must "as.character" to match the species name to tree labels in phyr
##otherwise the species names and trees don't match
dat1$sp = as.character(dat1$taxon)

##get phylogenies
phylo<-read.nexus("CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_resolvedUltra_July2020.nex")

##bayesian computation
##both species & island fixed effects

##full model, no imputation
##make formulae
##based on thse results get rid of problematic variables maxe & fhar
full<-as.formula(y ~ mass + mas2 + area + elev + volc + hurr + fore + loss + impa  + mong + (1|sp__) + (1|Island))

##use only one tree
modf<-pglmm(full, data = dat1, cov_ranef=list(sp = phylo[[50]]), family="binomial", bayes = T, marginal.summ="median", calc.DIC=T, prior= "pc.prior.auto")

##make a new dataset wth imputed values
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

##save the models
save.image("phyr100binomial.Rdata")

##significant full model
sigf<-as.formula(y ~ mass + mas2 + elev  + hurr + (1|sp__) + (1|Island))

##significant imputed model
sigi<-as.formula(y ~ mass + mas2 + elev + (1|sp__) + (1|Island))

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
save.image("phyr100binomial.Rdata")

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
save.image("phyr100binomial.Rdata")

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
save.image("phyr100binomial.Rdata")
