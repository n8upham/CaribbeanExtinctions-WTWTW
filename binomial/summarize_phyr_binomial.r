library(MASS)
library(geiger)
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
load("phyr100binomial.Rdata")

##missing data model
##get all mean sample wide coeffs into a single table
fix.m<-as.data.frame(sapply(save22, function(save22) save22$B,simplify=T))

##label rows, columns stay as is
rownames(fix.m)<-save22[[1]]$inla.model$names.fixed

##get extreme values
fix.ci<- as.data.frame(sapply(save22, function(save22) unlist(save22$B.ci),simplify=T))

##label rows, columns stay as is
rownames(fix.ci)<-paste(save22[[1]]$inla.model$names.fixed, rownames(fix.ci))

##put values together
fix<-rbind(fix.m, fix.ci)

##summarize sample wide coefficients
all.f<-as.data.frame(matrix(apply(fix, 1, median), ncol=3))

##label rows and cols
rownames(all.f)<-save22[[1]]$inla.model$names.fixed
colnames(all.f)<-c("mean", "0.025", "0.975")

##print out
sink("fixed_miss_data.txt")
print(all.f)
sink()

##add rownames as label
all.f$ID<-as.factor(rownames(all.f))

##reorder so they plot right
all.f$ID <- factor(all.f$ID, levels = c("hurr", "maxe", "elev" , "mas2", "mass", "(Intercept)"))

##create a new column with colors
all.f$col<-ifelse(all.f$'0.975'<0, "grey20", ifelse(all.f$'0.025'>0, "grey20", "grey80"))

##plot sample-wide effects
fix.plot<-ggplot(all.f,aes(ID, y=`mean`, ymin= `0.025`, ymax= `0.975`))+geom_hline(yintercept=0, linetype="dashed", colour="grey75", size=.75)+geom_pointrange(colour=all.f$col)+theme_pander()+ coord_flip()+ylab("Coefficient on probability of survival")+xlab("Sample-wide effects")+ scale_x_discrete(labels=c( "elev"="Mean island elevation", "hurr"="Hurricane index", "mass"="Mass", "mas2"= bquote(''*Mass^2*''), "(Intercept)"="Intercept"))
       
##save
ggsave("fixed_miss_data.pdf", h=5, w=5)

##non phylogenetic species effects
##summarize already by median
spe.m<-apply(as.data.frame(sapply(save22, function(save22) save22$inla.model$summary.random$`inla_effects[[1]]`[,c('mean')],simplify=T)), 1, median)

##non phylogenetic species effects extremes
spe.c1<-apply(as.data.frame(sapply(save22, function(save22) save22$inla.model$summary.random$`inla_effects[[1]]`[,c('0.025quant')],simplify=T)), 1, median)
spe.c2<-apply(as.data.frame(sapply(save22, function(save22) save22$inla.model$summary.random$`inla_effects[[1]]`[,c('0.975quant')],simplify=T)), 1, median)

##make into a single table
spe.r<-as.data.frame(cbind(spe.m, spe.c1, spe.c2))

##add rownames, ID and colnames
spe.r$ID<-save22[[1]]$inla.model$summary.random$`inla_effects[[1]]`[,c('ID')]
colnames(spe.r)<-colnames(all.f)[1:4]
rownames(spe.r)<-levels(save22[[1]]$random.effects$'1|sp'$sp)

##create a new column with colors
spe.r$col<-ifelse(spe.r$'0.975'<0, "grey20", ifelse(spe.r$'0.025'>0, "grey20", "grey80"))

##plot species effects
spe.plot<-ggplot(spe.r, aes(reorder(rownames(spe.r), desc(rownames(spe.r))), y=`mean`, ymin=`0.025`,ymax=`0.975`))+geom_hline(yintercept=0, linetype="dashed", colour="grey75", size=.75)+geom_pointrange(colour=spe.r$col)+theme_few()+ coord_flip()+ylab("Intercept on probability of survival")+xlab("Species")+ theme(axis.text.y = element_text(face="italic", size=5.5))
       
##save
ggsave("species_miss_data.pdf", h=6, w=6)

##island effects
##summarize already by median
isl.m<-apply(as.data.frame(sapply(save22, function(save22) save22$inla.model$summary.random$`inla_effects[[3]]`[,c('mean')],simplify=T)), 1, median)

##non phylogenetic species effects extremes
isl.c1<-apply(as.data.frame(sapply(save22, function(save22) save22$inla.model$summary.random$`inla_effects[[3]]`[,c('0.025quant')],simplify=T)), 1, median)
isl.c2<-apply(as.data.frame(sapply(save22, function(save22) save22$inla.model$summary.random$`inla_effects[[3]]`[,c('0.975quant')],simplify=T)), 1, median)

##make into a single table
isl.r<-as.data.frame(cbind(isl.m, isl.c1, isl.c2))

##add rownames, ID and colnames
isl.r$ID<-save22[[1]]$inla.model$summary.random$`inla_effects[[3]]`[,c('ID')]
colnames(isl.r)<-colnames(all.f)[1:4]
rownames(isl.r)<-levels(save22[[1]]$random.effects$'1|Island'$Island)

##create a new column with colors
isl.r$col<-ifelse(isl.r$'0.975'<0, "grey20", ifelse(isl.r$'0.025'>0, "grey20", "grey80"))

##plot island effects
isl.plot<-ggplot(isl.r, aes(reorder(rownames(isl.r), desc(rownames(isl.r))), y=`mean`, ymin=`0.025`,ymax=`0.975`))+geom_hline(yintercept=0, linetype="dashed", colour="grey75", size=.75)+geom_pointrange(colour=isl.r$col)+theme_few()+ coord_flip()+ylab("Intercept on probability of survival")+xlab("Islands")+ theme(axis.text.y = element_text(size=5))
       
##save
ggsave("islands_miss_data.pdf", h=8, w=4)

##imputed data model

##get all mean sample wide coeffs into a single table
fix.m<-as.data.frame(sapply(save24, function(save24) save24$B,simplify=T))

##label rows, columns stay as is
rownames(fix.m)<-save24[[1]]$inla.model$names.fixed

##get extreme values
fix.ci<- as.data.frame(sapply(save24, function(save24) unlist(save24$B.ci),simplify=T))

##label rows, columns stay as is
rownames(fix.ci)<-paste(save24[[1]]$inla.model$names.fixed, rownames(fix.ci))

##put values together
fix<-rbind(fix.m, fix.ci)

##summarize sample wide coefficients
all.f<-as.data.frame(matrix(apply(fix, 1, median), ncol=3))

##label rows and cols
rownames(all.f)<-save24[[1]]$inla.model$names.fixed
colnames(all.f)<-c("mean", "0.025", "0.975")

##print out
sink("fixed_impu_data.txt")
print(all.f)
sink()

##add rownames as label
all.f$ID<-as.factor(rownames(all.f))

##reorder so they plot right
all.f$ID <- factor(all.f$ID, levels = c("fhar", "elev" , "mas2", "mass", "(Intercept)"))

##create a new column with colors
all.f$col<-ifelse(all.f$'0.975'<0, "grey20", ifelse(all.f$'0.025'>0, "grey20", "grey80"))

##plot sample-wide effects
fix.plot<-ggplot(all.f,aes(ID, y=`mean`, ymin= `0.025`, ymax= `0.975`))+geom_hline(yintercept=0, linetype="dashed", colour="grey75", size=.75)+geom_pointrange(colour=all.f$col)+theme_pander()+ coord_flip()+ylab("Coefficient on probability of survival")+xlab("Sample-wide effects")+ scale_x_discrete(labels=c("mass"="Mass", "elev"="Mean island elevation", "mas2"= bquote(''*Mass^2*''), "(Intercept)"="Intercept"))
       
##save
ggsave("fixed_impu_data.pdf", h=5, w=5)

##non phylogenetic species effects
##summarize already by median
spe.m<-apply(as.data.frame(sapply(save24, function(save24) save24$inla.model$summary.random$`inla_effects[[1]]`[,c('mean')],simplify=T)), 1, median)

##non phylogenetic species effects extremes
spe.c1<-apply(as.data.frame(sapply(save24, function(save24) save24$inla.model$summary.random$`inla_effects[[1]]`[,c('0.025quant')],simplify=T)), 1, median)
spe.c2<-apply(as.data.frame(sapply(save24, function(save24) save24$inla.model$summary.random$`inla_effects[[1]]`[,c('0.975quant')],simplify=T)), 1, median)

##make into a single table
spe.r<-as.data.frame(cbind(spe.m, spe.c1, spe.c2))

##add rownames, ID and colnames
spe.r$ID<-save24[[1]]$inla.model$summary.random$`inla_effects[[1]]`[,c('ID')]
colnames(spe.r)<-colnames(all.f)[1:4]
rownames(spe.r)<-levels(save24[[1]]$random.effects$'1|sp'$sp)

##create a new column with colors
spe.r$col<-ifelse(spe.r$'0.975'<0, "grey20", ifelse(spe.r$'0.025'>0, "grey20", "grey80"))

##plot species effects
spe.plot<-ggplot(spe.r, aes(reorder(rownames(spe.r), desc(rownames(spe.r))), y=`mean`, ymin=`0.025`,ymax=`0.975`))+geom_hline(yintercept=0, linetype="dashed", colour="grey75", size=.75)+geom_pointrange(colour=spe.r$col)+theme_few()+ coord_flip()+ylab("Intercept on probability of survival")+xlab("Species")+ theme(axis.text.y = element_text(face="italic", size=5.5))
       
##save
ggsave("species_impu_data.pdf", h=6, w=6)

##island effects
##summarize already by median
isl.m<-apply(as.data.frame(sapply(save24, function(save24) save24$inla.model$summary.random$`inla_effects[[3]]`[,c('mean')],simplify=T)), 1, median)

##non phylogenetic species effects extremes
isl.c1<-apply(as.data.frame(sapply(save24, function(save24) save24$inla.model$summary.random$`inla_effects[[3]]`[,c('0.025quant')],simplify=T)), 1, median)
isl.c2<-apply(as.data.frame(sapply(save24, function(save24) save24$inla.model$summary.random$`inla_effects[[3]]`[,c('0.975quant')],simplify=T)), 1, median)

##make into a single table
isl.r<-as.data.frame(cbind(isl.m, isl.c1, isl.c2))

##add rownames, ID and colnames
isl.r$ID<-save24[[1]]$inla.model$summary.random$`inla_effects[[3]]`[,c('ID')]
colnames(isl.r)<-colnames(all.f)[1:4]
rownames(isl.r)<-levels(save24[[1]]$random.effects$'1|Island'$Island)

##create a new column with colors
isl.r$col<-ifelse(isl.r$'0.975'<0, "grey20", ifelse(isl.r$'0.025'>0, "grey20", "grey80"))

##plot island effects
isl.plot<-ggplot(isl.r, aes(reorder(rownames(isl.r), desc(rownames(isl.r))), y=`mean`, ymin=`0.025`,ymax=`0.975`))+geom_hline(yintercept=0, linetype="dashed", colour="grey75", size=.75)+geom_pointrange(colour=isl.r$col)+theme_few()+ coord_flip()+ylab("Coefficient on probability of survival")+xlab("Islands")+ theme(axis.text.y = element_text(size=5))
       
##save
ggsave("islands_impu_data.pdf", h=8, w=4)

##all predictors mdoel

##get all mean sample wide coeffs into a single table
fix.m<-as.data.frame(sapply(savema, function(savema) savema$B,simplify=T))

##label rows, columns stay as is
rownames(fix.m)<-savema[[1]]$inla.model$names.fixed

##get extreme values
fix.ci<- as.data.frame(sapply(savema, function(savema) unlist(savema$B.ci),simplify=T))

##label rows, columns stay as is
rownames(fix.ci)<-paste(savema[[1]]$inla.model$names.fixed, rownames(fix.ci))

##put values together
fix<-rbind(fix.m, fix.ci)

##summarize sample wide coefficients
all.f<-as.data.frame(matrix(apply(fix, 1, median), ncol=3))

##label rows and cols
rownames(all.f)<-savema[[1]]$inla.model$names.fixed
colnames(all.f)<-c("mean", "0.025", "0.975")

##print out
sink("fixed_full_data.txt")
print(all.f)
sink()

##add rownames as label
all.f$ID<-as.factor(rownames(all.f))

##reorder so they plot right
##fhar was excluded
all.f$ID <- factor(all.f$ID, levels = c("mong", "impa", "loss", "fore", "hurr", "volc", "elev" , "area", "mas2", "mass", "(Intercept)"))


##create a new column with colors
all.f$col<-ifelse(all.f$'0.975'<0, "grey20", ifelse(all.f$'0.025'>0, "grey20", "grey80"))

##plot sample-wide effects
fix.plot<-ggplot(all.f,aes(ID, y=`mean`, ymin= `0.025`, ymax= `0.975`))+geom_hline(yintercept=0, linetype="dashed", colour="grey75", size=.75)+geom_pointrange(colour=all.f$col)+theme_pander()+ coord_flip()+ylab("Coefficient on probability of survival")+xlab("Sample-wide effects")+ scale_x_discrete(labels=c("mong"="Mongoose present", "hurr"="Hurricane frequency", "volc"="Active volcano", "impa"="Human Footprint Index", "loss"="Forest loss (2000-2014)", "fore"="Forest cover (2000)", "elev"="Mean island elevation", "area"=" Island area", "mas2"= bquote(''*Mass^2*''), "mass"="Mass", "(Intercept)"="Intercept"))
       
##save
ggsave("fixed_full_data.pdf", h=6, w=6)

##non phylogenetic species effects
##summarize already by median
spe.m<-apply(as.data.frame(sapply(savema, function(savema) savema$inla.model$summary.random$`inla_effects[[1]]`[,c('mean')],simplify=T)), 1, median)

##non phylogenetic species effects extremes
spe.c1<-apply(as.data.frame(sapply(savema, function(savema) savema$inla.model$summary.random$`inla_effects[[1]]`[,c('0.025quant')],simplify=T)), 1, median)
spe.c2<-apply(as.data.frame(sapply(savema, function(savema) savema$inla.model$summary.random$`inla_effects[[1]]`[,c('0.975quant')],simplify=T)), 1, median)

##make into a single table
spe.r<-as.data.frame(cbind(spe.m, spe.c1, spe.c2))

##add rownames, ID and colnames
spe.r$ID<-savema[[1]]$inla.model$summary.random$`inla_effects[[1]]`[,c('ID')]
colnames(spe.r)<-colnames(all.f)[1:4]
rownames(spe.r)<-levels(savema[[1]]$random.effects$'1|sp'$sp)

##create a new column with colors
spe.r$col<-ifelse(spe.r$'0.975'<0, "grey20", ifelse(spe.r$'0.025'>0, "grey20", "grey80"))

##plot species effects
spe.plot<-ggplot(spe.r, aes(reorder(rownames(spe.r), desc(rownames(spe.r))), y=`mean`, ymin=`0.025`,ymax=`0.975`))+geom_hline(yintercept=0, linetype="dashed", colour="grey75", size=.75)+geom_pointrange(colour=spe.r$col)+theme_few()+ coord_flip()+ylab("Intercept on probability of survival")+xlab("Species")+ theme(axis.text.y = element_text(face="italic", size=5.5))
       
##save
ggsave("species_full_data.pdf", h=6, w=6)

##island effects
##summarize already by median
isl.m<-apply(as.data.frame(sapply(savema, function(savema) savema$inla.model$summary.random$`inla_effects[[3]]`[,c('mean')],simplify=T)), 1, median)

##non phylogenetic species effects extremes
isl.c1<-apply(as.data.frame(sapply(savema, function(savema) savema$inla.model$summary.random$`inla_effects[[3]]`[,c('0.025quant')],simplify=T)), 1, median)
isl.c2<-apply(as.data.frame(sapply(savema, function(savema) savema$inla.model$summary.random$`inla_effects[[3]]`[,c('0.975quant')],simplify=T)), 1, median)

##make into a single table
isl.r<-as.data.frame(cbind(isl.m, isl.c1, isl.c2))

##add rownames, ID and colnames
isl.r$ID<-savema[[1]]$inla.model$summary.random$`inla_effects[[3]]`[,c('ID')]
colnames(isl.r)<-colnames(all.f)[1:4]
rownames(isl.r)<-levels(savema[[1]]$random.effects$'1|Island'$Island)

##create a new column with colors
isl.r$col<-ifelse(isl.r$'0.975'<0, "grey20", ifelse(isl.r$'0.025'>0, "grey20", "grey80"))

##plot island effects
isl.plot<-ggplot(isl.r, aes(reorder(rownames(isl.r), desc(rownames(isl.r))), y=`mean`, ymin=`0.025`,ymax=`0.975`))+geom_hline(yintercept=0, linetype="dashed", colour="grey75", size=.75)+geom_pointrange(colour=isl.r$col)+theme_few()+ coord_flip()+ylab("Coefficient on probability of survival")+xlab("Islands")+ theme(axis.text.y = element_text(size=5))
       
##save
ggsave("islands_full_data.pdf", h=8, w=4)

##plot logistic regression result
## based on https://blogs.uoregon.edu/rclub/2016/04/05/plotting-your-logistic-regression-models/

##desired range of mass and mass^2
m_range <- seq(from=min(na.omit(dat1$mass)), to=max(na.omit(dat1$mass)), by=.01)

##X values 
impa_val <- mean(na.omit(dat1$impa))
loss_val <- mean(na.omit(dat1$loss))
fore_val <- mean(na.omit(dat1$fore))
hurr_val <- mean(na.omit(dat1$hurr))
elev_val <- mean(na.omit(dat1$elev))
area_val <- mean(na.omit(dat1$elev))

##make values
##reference group  volcano and mongoose set to 0
a_logits <- all.f['(Intercept)','mean'] + all.f['mass','mean'] * m_range + all.f['mas2','mean'] * m_range^2 + all.f['area','mean'] * area_val + all.f['elev','mean'] * elev_val + all.f['hurr','mean'] * hurr_val + all.f['fore','mean'] * fore_val + all.f['loss','mean'] * loss_val + all.f['impa','mean'] * impa_val

##mongoose but no volc
#b_logits <- all.f['(Intercept)','mean'] + all.f['mass','mean'] * m_range + all.f['mas2','mean'] * m_range^2 + all.f['area','mean'] * area_val + all.f['elev','mean'] * elev_val + all.f['hurr','mean'] * hurr_val + all.f['fore','mean'] * fore_val + all.f['loss','mean'] * loss_val + all.f['impa','mean'] * impa_val + all.f['mong','mean']

##volc but no mongoose
c_logits <- all.f['(Intercept)','mean'] + all.f['mass','mean'] * m_range + all.f['mas2','mean'] * m_range^2 + all.f['area','mean'] * area_val + all.f['elev','mean'] * elev_val + all.f['hurr','mean'] * hurr_val + all.f['fore','mean'] * fore_val + all.f['loss','mean'] * loss_val + all.f['impa','mean'] * impa_val + all.f['volc','mean']

# Compute the probibilities (this is what will actually get plotted):
a_probs <- exp(a_logits)/(1 + exp(a_logits))
#b_probs <- exp(b_logits)/(1 + exp(b_logits))
c_probs <- exp(c_logits)/(1 + exp(c_logits))

# first you have to get the information into a long dataframe, which is what ggplot likes :)
plot.data <- data.frame(baseline=a_probs, volcano=c_probs, mast=m_range)
plot.data <- gather(plot.data, key=group, value=prob, baseline:volcano)

##print out the table to get probabilities
write.csv(plot.data, "all_predict_table.csv")

##transform mass back
plot.data$mass<-10^(plot.data$mast*sd(log10(na.omit(dat1$Body_mass))) + mean(log10(na.omit(dat1$Body_mass))))

##actual observations of body mass
x<-dat1$Body_mass[!is.na(dat1$Body_mass)]

##plot in original values but showing log scale
m.plot<-ggplot(plot.data, aes(x=mass, y=prob, color=group)) + # asking it to set the color by the variable "group" is what makes it draw three different lines
  geom_line(lwd=1.5) + 
  labs(x="Mass in Kg", y="Probability of survival")+theme_few()+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+scale_colour_grey()+geom_rug(data =as.data.frame(x), inherit.aes = F, aes(x=x))+ theme(legend.title = element_blank(), legend.position="top")

##save
ggsave("mass_full_data.pdf", h=4, w=6)

##remove all
rm(list=ls()) 
