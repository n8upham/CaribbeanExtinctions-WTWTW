library(MASS)
library(INLA)
library(parallel)
library(phyr)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
library(scales)

##summarizing inla results across 100 trees
##remove prior data
rm(list=ls()) 

##get the models resulting from phylo_trees_phyr.r
load("phyr100binomial_v3.Rdata")

##get imputed values for mass
mass.m<-modf$inla.model$summary.fitted.values[mass.na, ]

##scale back to kilograms
mass.m<-10^(mass.m*sd(log10(na.omit(dat1$Body_mass))) + mean(log10(na.omit(dat1$Body_mass))))

##print
sink("Hexolobodontinae_x_v2.txt")
print("mass as entered")
print(modf$inla.model$summary.fitted.values[mass.na, ])
print("mass in Kg")
print(mass.m)
sink()

##remove all
rm(list=ls()) 
