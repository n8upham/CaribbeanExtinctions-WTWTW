library(geiger)
library(phyr)
library(INLA)
library(ppcor)
library(reshape2)
library(corrplot)
library("psych")

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

##test for multicollinearity
##missing value rows excluded in pairs as os opposed to all rows

## focus only on island data characteristics
dat3<-dat1[!duplicated(dat1[,'Island']),]

##get corr coefficients 
mcor<- cor(dat3[,27:36], method = "pearson", use="pairwise.complete.obs")
diag(mcor)<-NA	
mcor.m<-na.omit(melt(mcor, value.name="corr"))

##get rid of duplicates
mcor.m<-mcor.m[!duplicated(mcor.m[,'corr']),]

##calculate pairwise correlation p-values using an inelegant loop
p<-list()

for(i in 1:dim(mcor.m)[1]){p[[i]]<-cor.test(dat3[, as.character(mcor.m$Var1[i])], dat3[, as.character(mcor.m$Var2[i])], method="pearson")$p.value}

##add values to matrix
mcor.m$pval<-unlist(p)

##get vars highly and significantly correlated
m1.m<-subset(mcor.m, pval<0.05 & abs(corr) >0.7)

##order by corr
m1.m<-m1.m[order(m1.m$corr),]

##export for decision making
write.csv(m1.m, "collinear.csv", row.names =F)

##plot for visualization
colnames(mcor)<-c("area", "elevation", "maximum elevation", "volcano", "hurricane", "forest cover", "recent forest loss", "human impact", "first human arrival", "mongoose") 
rownames(mcor)<-colnames(mcor)

##matrix of p values
pcor<-corr.test(dat3[,27:36], method = "pearson", use="pairwise.complete.obs")$p
diag(pcor)<-NA

colnames(pcor)<-c("area", "elevation", "maximum elevation", "volcano", "hurricane", "forest cover", "recent forest loss", "human impact", "first human arrival", "mongoose") 
rownames(pcor)<-colnames(pcor)
	
pdf(file="Island_correlation.pdf")
corrplot(mcor, type = "lower", order = "original", tl.col = "black", tl.srt = 45, addCoef.col=TRUE, p.mat=pcor, sig.level = 0.01, diag=F)
dev.off()

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
full<-as.formula(y ~ mass + mas2 + area + elev + volc + hurr + fhar + fore + loss + impa  + mass:fhar + (1|sp__) + (1|Island))

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

##save the models
save.image("phyr100binomial_dat.Rdata")

