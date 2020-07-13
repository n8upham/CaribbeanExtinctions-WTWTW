library(MASS)
library(geiger)
library(MCMCglmm)
library(MCMCvis)
library(brms)
library(parallel)
library(mice)

##remove prior data
rm(list=ls()) 

##taxa already matching perfectly using treedata
##get phylogenies
tree<-read.nexus("CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_resolvedUltra_July2020.nex")

##detect cores
numCores <- detectCores()

##speed through making ultrametric
mkul <- function(i) {
  chronos(tree[[i]], lambda=0.5)
}

system.time({
  library(ape)
  phylo <- mclapply(1:100, mkul, mc.cores =numCores)
})

##rm attributes
for (i in 1:length(phylo)){attr(phylo[[i]],"order") <- NULL}
for (i in 1:length(phylo)){attr(phylo[[i]],"call") <- NULL}
for (i in 1:length(phylo)){attr(phylo[[i]],"ploglik") <- NULL}
for (i in 1:length(phylo)){attr(phylo[[i]],"rates") <- NULL}
for (i in 1:length(phylo)){attr(phylo[[i]],"message") <- NULL}
for (i in 1:length(phylo)){attr(phylo[[i]],"PHIIC") <- NULL}
for (i in 1:length(phylo)){attr(phylo[[i]],"class") <- "phylo"}
class(phylo) <- "multiPhylo"	#specify class

##save so these can be used
save.image("ultra_trees.Rdata")
