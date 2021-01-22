#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Code - Caribbean mammals from MamPhy 
###
# Movies - animated gifs showing tree uncertainty in the credible sets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####
# LOAD the 100 mammal trees of interest.
###
dirname<-"/Users/nate/Desktop/Carib_Extinctions/Turvey_Davalos_correlatesOfExtinction/CaribbeanExtinctions-WTWTW/mamPhy_pruningCode"
setwd(dirname)
library(ape); library(phyloch); library(phytools)

mamPhy_100<-read.nexus(file="CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_resolvedUltra_FINAL.nex")

# get root heights
subRoots<-c()
for(i in 1:length(mamPhy_100)){
	subRoots[i]<-max(node.depth.edgelength(phy=mamPhy_100[[i]]))
}
maxOfAll<-max(subRoots)
maxOfAll_ID<-match(maxOfAll,subRoots)

# which are extinct?
extinctSp<-na.omit(read.delim("nameToAdd_mrca1_mrca2_MamPhy_noParen.txt"))$nameToAdd

# set node colors
library(viridis); library(gplots); #library(diversitree)

RESDIR<-"gif_100trees"

# NOW ANIMATE
install.packages("magick")
library(magick)

img <- image_graph(height=800, width=600, res = 96, bg="white")
#pdf(file=paste0(RESDIR,"/credibleSet_CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_resolvedUltra_FINAL.pdf"),height=8,width=6, onefile = TRUE)

for(i in 1:100){

	tree_i<-ladderize(mamPhy_100[[i]])

	#plot oldest tree and axis
	plot(ladderize(mamPhy_100[[maxOfAll_ID]]), edge.width=0.5, show.tip.label=FALSE, cex=0.3, edge.color=NA, type="phylogram")
	axis(side=1,cex.lab=0.4,at=(maxOfAll-c(0,20,40,60,80)),labels=(c(0,20,40,60,80)))

	subRoot_BASE<-max(node.depth.edgelength(mamPhy_100[[maxOfAll_ID]]))

	rootTime<-100

	#plot target tree
	## get last phylo plot parameters
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	subRoot_TARGET<-max(node.depth.edgelength(tree_i))
	par(new=TRUE)
	plot((tree_i), edge.width=0.5, show.tip.label=FALSE, cex=0.3, edge.color="black", type="phylogram", 
		x.lim=c(subRoot_TARGET-obj$x.lim[2],subRoot_TARGET) )

	## get last phylo plot parameters
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	x.tip<-obj$xx[1:obj$Ntip]
	y.tip<-obj$yy[1:obj$Ntip]
	#tipsExtinct<-as.vector(na.omit(match(extinctSp,tree_i$tip.label)))#c(1:5911)[tree_i$tip.label %in% missingSp]
	#tipsExtinct<-y.tip[tree_i$tip.label %in% extinctSp]
	tipsExtinct<-c(1:obj$Ntip)[tree_i$tip.label %in% extinctSp]

	#label the TIPS for MISSING
	#tiplabels(cex=0.1)
	tiplabels(tip=tipsExtinct, frame="none",pch=23, bg="yellow",col="black", cex=1, offset=2) # 
	#lines(y=cbind(tipsExtinct,tipsExtinct), x=cbind(rep(x.tip[1],length(tipsExtinct)), rep(x.tip[1]+5,length(tipsExtinct)) ), col="red", lty=1) # 
	#lines(y=tipsExtinct, x=rep(x.tip[1],length(tipsExtinct)), col="red", lty=1) # 

	mtext(side=3, text=paste0("tree ",i), font=2, cex=0.9)
	mtext(side=1, text="Millions of years before present (Ma)", cex=1, line=3)

}
dev.off()

animation <- image_animate(img, fps = 10)
print(animation)
image_write(animation, paste0(RESDIR,"/credibleSet_credibleSet_CaribbeanMam_77taxa_MamPhy_added42taxa_random100trees_resolvedUltra_FINAL.gif"))



