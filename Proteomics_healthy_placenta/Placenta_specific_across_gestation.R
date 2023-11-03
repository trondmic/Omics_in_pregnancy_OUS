###############
#
#
# Author: Maren-Helene Langeland Degnes
# Purpose: Define placenta specific proteins and
#           model these across gestation
#
# Lines 109 to 121 in this script is based on an original R script developed 
# by A. Tarca (atarca@med.wayne.edu), provided to us by the author on 22nd of November 2020
# Reference: 
#   Erez, O., Romero, R., Maymon, E., Chaemsaithong, P., 
#   Done, B., Pacora, P., Panaitescu, B., Chaiworapongsa, 
#   T., Hassan, S. S., & Tarca, A. L. (2017). The prediction 
#   of late-onset preeclampsia: Results from a longitudinal 
#   proteomics study. PLOS ONE, 12(7), e0181468. 
#
#############





# Load packages
library(tidyverse)
library(Hotelling)
library(lme4)
library(splines)
library(dendextend)
library(ape)
library(RColorBrewer)
library(cluster)
library(scales)
library(openxlsx)





################             Load data            ###################
# load 4-vessel data set median normalized venoarterial differences 
load("data/MTT_FBmedians_placentanclock.RData")
# load STORK data set 
load("data/Oslo_STORK_log2win_exprset.RData")
samprot <- exprs(Oslo_STORK_log2win_exprset) %>% t()
samp <- pData(Oslo_STORK_log2win_exprset)
prot <- fData(Oslo_STORK_log2win_exprset)








######################### Define placenta specific proteins ##################
# upregulated proteins are
#             1. higher secreted from placenta to mother than in arm
#
upr.idx <- which(MatNormTest$`adjusted p`<0.05
                 & MatNormTest$`t statistic`>0
                 & MatNormTest$PLARMpadj<0.05
                 & MatNormTest$PLARMt>0)

# specific proteins are
#             1. secreted from placenta to mother
#             2. AND difference in arm should be 0 or negative
spec.idx <- which(MatNormTest$`adjusted p`<0.05
                  & MatNormTest$`t statistic`>0
                  & ((ArmNormTest$`p value`>0.05)
                     | ArmNormTest$`t statistic`<0))

# Save the three different categories
spec_but_not_upr <- setdiff(spec.idx,upr.idx)
spec_and_upr <- intersect(spec.idx,upr.idx)
upr_but_not_spec <- setdiff(upr.idx,spec.idx)
prot.indeces <- c(spec_but_not_upr,spec_and_upr,upr_but_not_spec) %>% unique()


# Save information about selected proteins for short access later
prodclass.gen <- ifelse(1:nrow(MatNormTest)%in%spec_but_not_upr,1,ifelse(1:nrow(MatNormTest)%in%upr_but_not_spec,3,2))
prodclass <- prodclass.gen[prot.indeces]
prodclass.protnames <- prot$TargetFullName[prot.indeces]








###########################  Mixed model per protein ##########################
# Create objects 
ST.pred.mx <- c()
pred.mx <- c()
models <- list()
modelsi <- 1
pdfpfx <- "UprSpec_dev"
bsreslist <- list()   # container for testing coefficients in mm

# Save objects 
normal_longit.idx <- which(!is.na(samp$TimePoint) & samp$Control==1)
ST.ID <- samp$ID[normal_longit.idx]                          # ID for only stork cohort
meanTP <- tapply(samp$GAWeeks[normal_longit.idx] %>% as.numeric(),samp$TimePoint[normal_longit.idx],mean)

# Model and plot all proteins indexed inside prot.indeces
for (prnum in 1:length(prot.indeces)){
  par(mar=c(3,3,2,0.5),xpd=T,cex.main=0.9,ps=8)
  pr <- prot.indeces[prnum]
  # set data
  ST.proteinmeas <- samprot[c(normal_longit.idx),pr]
  proteinname <- prot$TargetFullName[pr]
  uniprot <- prot$UniProt[pr]
  aptname <- prot$AptName[pr]
  
  
  ####  mixed modelling for only stork samples
  x <- samp$GAWeeks[normal_longit.idx] %>% as.numeric()
  y <- ST.proteinmeas
  bsplines <- bs(x, degree=3,knots=meanTP[c(2)], Boundary.knots =meanTP[c(1,3)])
  bsplines.df <- bsplines %>% as.data.frame()
  colnames(bsplines.df) <- paste0("bs",1:4)
  mmfit <- lmer(y ~ 1 + bs1 + bs2 + bs3 + bs4 + (1| ST.ID),
                data=cbind(bsplines.df,ST.ID,y),REML=F)            #x = ST.GAWeeks, y = observed proteinmeas
  bsreslist[[prnum]] <- lmerTest::as_lmerModLmerTest(mmfit) %>% summary()
  xnew <- seq(meanTP[1],meanTP[3],length.out = 100)
  xnew.df <- predict(bsplines, xnew) %>% as.data.frame()
  colnames(xnew.df) <- paste0("bs",1:4)
  newdata <- cbind(xnew.df)
  ynew <- predict(mmfit, newdata=newdata, re.form=~0)                     #ynew = predicted proteinmeas at xnew
  ST.pred.mx <- rbind(ST.pred.mx,ynew)
  
}






##################### Cluster each protein's estimated mean curve ######################
# Create dendrogram objects
cs <- 1-(ST.pred.mx %>% t() %>% cor(.,method="pearson"))
csd <- cs %>% as.dist(); hc <- csd %>% hclust(); hc$labels <- paste0(prot$TargetFullName[prot.indeces],"-",prot$EntrezGeneSymbol[prot.indeces])
dend <- hc %>% as.dendrogram()

# Detect optimal number of clusters above two clusters
silh.means <- c()
for(ncl in 2:10){silh.means <- c(silh.means,summary(silhouette(cutree(hc,k=ncl),csd))$avg.width)}
par(mfrow=c(1,1),mar=(rep(4,4)));plot(2:10,silh.means,type="b",
                                      main="Mean silhouette widths for each number of clusters 2:10",
                                      ylab="Mean silhouette widths")
optimal.num.cl <-  grep(max(silh.means),silh.means) + 1

# Cut tree by a cutoff and save colors for plotting
clusMember <- cutree(dend,optimal.num.cl)
allcolors <- c(brewer.pal(n=8,name="Dark2"),"darkcyan","firebrick1")
clusColors <- allcolors[1:optimal.num.cl]
labels_colors(dend) <- clusColors[clusMember][order.dendrogram(dend)]










################==================== Combined plot of dendrogram and scaled developments in pdf ==================###############
pdf('results/Figure 3.pdf',
    width=6.6, height=4.5)

# Set all sizes and the plot window layout
dotcex <- 0.5
leavescex <- 0.5
axisvaluesdendcex <- 1
axislabelscex <- 0.5
plotmaincex <- 1
legendcex <- 0.5
linethickness <- 0.5
xlabelcex <- 0.7
layout(mat=matrix(c(rep(1,optimal.num.cl*2),2:(optimal.num.cl+1)),nrow=3,byrow=TRUE))
par(mar=c(15,3,0.4,0),
    ps=8,
    cex.axis=axisvaluesdendcex,
    cex.lab=axislabelscex,
    cex.main=plotmaincex,
    mgp=c(3,0.5,0.5))

# Adjust dendrogram settings
dend <- dend %>% 
  set("labels_cex",0.7) %>% 
  set("branches_lwd", 0.5)

# Create cluster name vector
clusternames <- paste("Cluster", c(2,1)) #c(2,6,5,1,4,3) #c(6, 1, 3, 4, 2,5)

# plot dendrogram
plot(dend,edgePar = list(t.col=clusColors),ylab="",yaxt="n")
mtext("Pearson correlation",side=2,line=2,outer=F,cex=xlabelcex,at=1)
axis(2,at=seq(0,2,by=0.5),labels=seq(1,-1,by=-0.5))
abline(h=1.75,lty=2,lwd=0.5)
text(-2,2.01,"A.",cex=1,xpd=NA,face="bold")

# Settings for development curves beneath dendrogram
par(mar=c(1,1,1,1))
clusterdata <- ST.pred.mx
x <- seq(meanTP[1],meanTP[3],length.out = 100)
firstplot <- T
for (cluster in c(2,1) 
){
  clustercol <- clusColors[cluster]
  clus <- which(clusMember==cluster)
  
  clusterdata[clus,] <- clusterdata[clus,] + (mean(clusterdata[clus,1]) - clusterdata[clus,1])
  clusterdata[clus,] <- clusterdata[clus,] - mean(clusterdata[clus,1])
  
  if(!firstplot){
    if(cluster==4){mtext("--Gestational age in weeks--",side=1,line=2,outer=F,at=23,cex=xlabelcex)}
    par(mar=c(3,2,1,1))
    plot(x,clusterdata[clus[1],],type="l",
         ylim=c(min(clusterdata[clus,]),max(clusterdata[clus,])),col=clustercol,
         xlab="",ylab="",
         lwd=linethickness,
         main=clusternames[cluster],
         xlim=c(15,32))
  }else{
    par(mar=c(3,3,1,1))
    plot(x,clusterdata[clus[1],],type="l",
         ylim=c(min(clusterdata[clus,]),max(clusterdata[clus,])),col=clustercol,
         xlab="",ylab="",lwd=linethickness,
         main=clusternames[cluster],
         xlim=c(15,32))
    mtext("Scaled log2 RFU",side=2,line=2,outer=F,cex=xlabelcex)
    
    text(15.3,2.8,"B.",cex=1, xpd=NA,ps=5)
  }
  if (length(clus)>1){
    for (protein in clus[2:length(clus)]){
      lines(x,clusterdata[protein,],col=clustercol,lwd=linethickness,xlim=c(15,32))
    }
  }
  firstplot <- F
}




dev.off()


