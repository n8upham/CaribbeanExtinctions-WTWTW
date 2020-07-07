###########
# Code to prune -and impute- MAMMAL trees to Caribbean taxa

# Load packages:
library(ape); library(phytools); library(PDcalc); library(phangorn)
#install.packages("remotes")
#remotes::install_github("davidnipperess/PDcalc")

#****************
# Upham et al. 2019 trees
#####
setwd("/Users/Nate/Desktop/Carib_Extinctions/Turvey_Davalos_correlatesOfExtinction/mamPhy_pruningCode_final")
trees100<-read.nexus(file="MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_v2_sample100_nexus.trees")

taxonList<-na.omit(read.table("MamPhy_CaribbeanTaxaToPruneOut.txt", header=TRUE))

toDrop<-setdiff(trees100[[1]]$tip.label, as.vector(taxonList$MamPhy_tip))
toRename<-as.vector(taxonList$nameToAdd)

trees100_carib<-vector("list",length(trees100))
for (i in 1:length(trees100)){
    dropped<-drop.tip(trees100[[i]],toDrop)
    dropped$tip.label<-toRename[match(dropped$tip.label,as.vector(taxonList$MamPhy_tip))]
    trees100_carib[[i]]<-dropped
}
write.nexus(trees100_carib,file="CaribbeanMam_35taxa-plus3_MamPhy_random100trees.nex")

pdf(file="plottedTree_CaribbeanMam_35taxa-plus3_MamPhy_1of100.pdf")
plot(trees100_carib[[1]],cex=0.6, label.offset=1)
#plot(drop.tip(trees100[[1]],toDrop),cex=0.6, label.offset=1)

nodelabels(cex=0.5)
axisPhylo()
dev.off()

taxaToAdd<-na.omit(read.table("nameToAdd_mrca1_mrca2_MamPhy.txt", header=TRUE))

newCaribTree100<-vector("list",length(trees100_carib))
for (i in 1:length(trees100_carib)){
caribTree<-trees100_carib[[i]]

    for (j in 1:length(taxaToAdd[,1])){
        toAdd<-as.vector(taxaToAdd$nameToAdd[j])
        mrca1<-as.vector(taxaToAdd$mrca1[j])
        mrca2<-as.vector(taxaToAdd$mrca2[j])
        node<-getMRCA(phy=caribTree,tip=c(mrca1, mrca2))
        depth<-max(nodeHeights(caribTree))-nodeheight(caribTree,node) # attaches the new tips to non-ultrametric trees at the same distance as the current tips at maximum distance from the root.

        caribTree<-bind.tip(caribTree,tip.label=toAdd,edge.length=depth,where=node,position=0)
    }
    # prune out the 2 extra sloths, and 1 extra primate
    caribTree<-drop.tip(caribTree, c("Choloepus_hoffmanni", "Choloepus_didactylus", "Callicebus_cupreus"))
    newCaribTree100[[i]]<-caribTree
}
write.nexus(newCaribTree100,file="CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees.nex")

pdf(file="plottedTree_CaribbeanMam_77taxa_MamPhy_added42taxa_all100.pdf", onefile=TRUE)
for(i in 1:length(newCaribTree100)){
    plot(ladderize(newCaribTree100[[i]]),cex=0.4, label.offset=1)
    nodelabels(cex=0.4)
    axisPhylo()
    mtext(side=3,text=paste("tree ",i,sep=""))
}
dev.off()

##
# Now RESOLVE those polytomies
resolvedCaribTree100<-vector("list",length(trees100_carib))
for (i in 1:length(trees100_carib)){
    caribTree<-newCaribTree100[[i]]
    resolvedTree<-multi2di(caribTree,random=TRUE)
        # change zero branch lengths to 0.0001
        resolvedTree$edge.length[resolvedTree$edge.length == 0] <- 0.0001

    # and MAKE ULTRAMETRIC:    
    readyTree<-nnls.tree(cophenetic(resolvedTree),resolvedTree,rooted=TRUE)
    is.ultrametric(readyTree) ## should pass

    resolvedCaribTree100[[i]]<-readyTree
}
write.nexus(resolvedCaribTree100,file="CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_resolvedUltra_FINAL.nex")

pdf(file="plottedTree_CaribbeanMam_77taxa_MamPhy_added42taxa_all100_resolvedUltra_FINAL.pdf", onefile=TRUE)
for(i in 1:length(resolvedCaribTree100)){
    plot(ladderize(resolvedCaribTree100[[i]]),cex=0.4, label.offset=1)
    #nodelabels(cex=0.4)
    axisPhylo()
    mtext(side=3,text=paste("tree ",i,sep=""))
}
dev.off()



