####### script information ##########
#
#  Author: Maren-Helene Langeland Degnes  
#  Date: 22.Aug.2022
#  Purpose: 
#           A. median normalization of venoarterial differences in placenta
#           B. winsorization of the median normalized venoarterial differences
#
#######

#### load packages ########
library(tidyverse)   # for piping
library(haven)       # to read sav-files (spss-files)
library(pls)     
library(scales)
library(openxlsx)
library(limma)













#################### A. Median normalization ###################

################## ============= Calculate median normalized values in control group ===================== ############################
rm(list=ls())
load(file="dataset.RData") # we get samp = phenodata, samprot = expression data, and prot = featuredata

#check that objects with pehotype data and expression data (Biobase nomenclature) match
(rownames(samp) == rownames(samprot)) %>% table()

## subset into Controls
subset.idx <- which(samp$PE==0)
samp <- samp[subset.idx,]
samprot <- samprot[subset.idx,]

############## ================ Choose Control patients from 4V and create new variables ============ ############
# subset into specific vessels
UE.samp <- samp[which(samp$SampleGroup=="UE"),]
UE.samprot <- samprot[which(samp$SampleGroup=="UE"),]

PAE.samp <- samp[which(samp$SampleGroup=="PAE"),]
PAE.samprot <- samprot[which(samp$SampleGroup=="PAE"),]

PVE.samp <- samp[which(samp$SampleGroup=="PVE"),]
PVE.samprot <- samprot[which(samp$SampleGroup=="PVE"),]

AUE.samp <- samp[which(samp$SampleGroup=="AUE"),]
AUE.samprot <- samprot[which(samp$SampleGroup=="AUE"),]

VUE.samp <- samp[which(samp$SampleGroup=="VUE"),]
VUE.samprot <- samprot[which(samp$SampleGroup=="VUE"),]

######### ========= Check that patients have the same order so they can be compared and normalized pairwise ======== #######
(UE.samp$ID == PAE.samp$ID) %>% table()
(VUE.samp$ID == AUE.samp$ID) %>% table()
(PVE.samp$ID == PAE.samp$ID) %>% table() 

############ =============== Calculate differences between pairs of vessels and subtract median ============= ###########
m.c.diffs <- UE.samprot-PAE.samprot
adj.m.c.diffs <- m.c.diffs
for (patient in 1:nrow(m.c.diffs)){ 
  adj.m.c.diffs[patient,] <- m.c.diffs[patient,] - median(m.c.diffs[patient,]) 
}

f.c.diffs <- VUE.samprot-AUE.samprot
adj.f.c.diffs <- f.c.diffs
for (patient in 1:nrow(f.c.diffs)){ 
  adj.f.c.diffs[patient,] <- f.c.diffs[patient,] - median(f.c.diffs[patient,]) 
}


marm.c.diffs <- PVE.samprot[which(PVE.samp$ID  %in% PVEPAE.patients),] - PAE.samprot[which(PAE.samp$ID  %in% PVEPAE.patients),]
adj.marm.c.diffs <- marm.c.diffs
for (patient in 1:nrow(marm.c.diffs)){ 
  adj.marm.c.diffs[patient,] <- marm.c.diffs[patient,] - median(marm.c.diffs[patient,]) 
}

############# ============== control! ============== ################
(prot$AptName == colnames(samprot)) %>% table()
(prot$AptName == colnames(UE.samprot)) %>% table()
(prot$AptName == colnames(adj.m.c.diffs)) %>% table()
(prot$AptName == colnames(m.c.diffs)) %>% table()
(prot$AptName == colnames(adj.marm.c.diffs)) %>% table()
(prot$AptName == colnames(marm.c.diffs)) %>% table()
(prot$AptName == colnames(adj.f.c.diffs)) %>% table()
(prot$AptName == colnames(f.c.diffs)) %>% table()
(rownames(samp) == rownames(samprot)) %>% table()
(apply(adj.m.c.diffs,1,median) == 0) %>% table()
(apply(adj.marm.c.diffs,1,median) == 0) %>% table()
(apply(adj.f.c.diffs,1,median) == 0) %>% table()

################# ============= Calculate median normalized values in PE group 
load(file="data/lgRFU.RData")

(rownames(samp) == rownames(samprot)) %>% table()

## subset into PE
subset.idx <- which(samp$PE==1)
samp <- samp[subset.idx,]
samprot <- samprot[subset.idx,]


############## ================ Choose Control patients from 4V and create new variables 
UE.samp <- samp[which(samp$SampleGroup=="UE"),]
UE.samprot <- samprot[which(samp$SampleGroup=="UE"),]

PAE.samp <- samp[which(samp$SampleGroup=="PAE"),]
PAE.samprot <- samprot[which(samp$SampleGroup=="PAE"),]

PVE.samp <- samp[which(samp$SampleGroup=="PVE"),]
PVE.samprot <- samprot[which(samp$SampleGroup=="PVE"),]

AUE.samp <- samp[which(samp$SampleGroup=="AUE"),]
AUE.samprot <- samprot[which(samp$SampleGroup=="AUE"),]

VUE.samp <- samp[which(samp$SampleGroup=="VUE"),]
VUE.samprot <- samprot[which(samp$SampleGroup=="VUE"),]

######### ========= Check that patients have the same order so they can be compared and normalized pairwise
(UE.samp$ID == PAE.samp$ID) %>% table()
(VUE.samp$ID == AUE.samp$ID) %>% table()
(PVE.samp$ID == PAE.samp$ID) %>% table() 

############ =============== Calculate differences between pairs of vessels and subtract median 
m.pe.diffs <- UE.samprot-PAE.samprot
adj.m.pe.diffs <- m.pe.diffs
for (patient in 1:nrow(m.pe.diffs)){ 
  adj.m.pe.diffs[patient,] <- m.pe.diffs[patient,] - median(m.pe.diffs[patient,]) 
}

f.pe.diffs <- VUE.samprot-AUE.samprot
adj.f.pe.diffs <- f.pe.diffs
for (patient in 1:nrow(f.pe.diffs)){ 
  adj.f.pe.diffs[patient,] <- f.pe.diffs[patient,] - median(f.pe.diffs[patient,]) 
}


marm.pe.diffs <- PVE.samprot[which(PVE.samp$ID  %in% PVEPAE.patients),] - PAE.samprot[which(PAE.samp$ID  %in% PVEPAE.patients),]
adj.marm.pe.diffs <- marm.pe.diffs
for (patient in 1:nrow(marm.pe.diffs)){ 
  adj.marm.pe.diffs[patient,] <- marm.pe.diffs[patient,] - median(marm.pe.diffs[patient,]) 
}


############# ============== control! 
(prot$AptName == colnames(samprot)) %>% table()
(prot$AptName == colnames(UE.samprot)) %>% table()
(prot$AptName == colnames(adj.m.pe.diffs)) %>% table()
(prot$AptName == colnames(m.pe.diffs)) %>% table()
(prot$AptName == colnames(adj.marm.pe.diffs)) %>% table()
(prot$AptName == colnames(marm.pe.diffs)) %>% table()
(prot$AptName == colnames(adj.f.pe.diffs)) %>% table()
(prot$AptName == colnames(f.pe.diffs)) %>% table()
(rownames(samp) == rownames(samprot)) %>% table()
(apply(adj.m.pe.diffs,1,median) == 0) %>% table()
(apply(adj.marm.pe.diffs,1,median) == 0) %>% table()
(apply(adj.f.pe.diffs,1,median) == 0) %>% table()




##############  Write r data file of datasets
load(file="data/lgRFU.RData")
save(samprot,
     samp,
     prot,
     adj.m.c.diffs,
     adj.marm.c.diffs,
     adj.f.c.diffs,
     adj.m.pe.diffs,
     adj.marm.pe.diffs,
     adj.f.pe.diffs,
     file="data/lgRFUmedians.RData")
















################### B. Winsorization of 4V differences #####################



###### winsorize differences in both tales
rm(list=ls())
load(file="data/lgRFUmedians.RData")

matdiffs <- rbind(adj.m.c.diffs,adj.m.pe.diffs)
fetdiffs <- rbind(adj.f.c.diffs,adj.f.pe.diffs)
armdiffs <- rbind(adj.marm.c.diffs,adj.marm.pe.diffs)

matdiffs <- apply(matdiffs,2,function(x){q<-quantile(x,c(0.99,0.01));uoutls<-which(x>q[1]);loutls<-which(x<q[2]); x[uoutls]<-q[1];x[loutls]<-q[2]; x })
fetdiffs <- apply(fetdiffs,2,function(x){q<-quantile(x,c(0.99,0.01));uoutls<-which(x>q[1]);loutls<-which(x<q[2]); x[uoutls]<-q[1];x[loutls]<-q[2]; x })
armdiffs <- apply(armdiffs,2,function(x){q<-quantile(x,c(0.99,0.01));uoutls<-which(x>q[1]);loutls<-which(x<q[2]); x[uoutls]<-q[1];x[loutls]<-q[2]; x })

mat.group.vec <- c(rep("C",nrow(adj.m.c.diffs)),rep("PE",nrow(adj.m.pe.diffs)))
adj.m.c.diffs <- matdiffs[which(mat.group.vec=="C"),]
adj.m.pe.diffs <- matdiffs[which(mat.group.vec=="PE"),]

fet.group.vec <- c(rep("C",nrow(adj.f.c.diffs)),rep("PE",nrow(adj.f.pe.diffs)))
adj.f.c.diffs <- fetdiffs[which(fet.group.vec=="C"),]
adj.f.pe.diffs <- fetdiffs[which(fet.group.vec=="PE"),]

marm.group.vec <- c(rep("C",nrow(adj.marm.c.diffs)),rep("PE",nrow(adj.marm.pe.diffs)))
adj.marm.c.diffs <- armdiffs[which(marm.group.vec=="C"),]
adj.marm.pe.diffs <- armdiffs[which(marm.group.vec=="PE"),]


save(samprot,
     samp,
     prot,
     adj.m.c.diffs,
     adj.marm.c.diffs,
     adj.f.c.diffs,
     adj.m.pe.diffs,
     adj.marm.pe.diffs,
     adj.f.pe.diffs,
     file="data/lgRFUmedianswins.RData")



sessionInfo()
