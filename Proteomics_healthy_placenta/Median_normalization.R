###############
#
#
# Author: Maren-Helene Langeland Degnes
# Purpose: Median normalize venoarterial differences
#
#
#############




# Load packages
rm(list=ls())
library(tidyverse)
library(haven)
library(pls)


# Load data - the 4-vessel data set with each vessel
f.path.to.dta <- "data/lgRFU_placentalclock.RData"
load(file=f.path.to.dta)







############## ================ Choose participants from 4V and create new variables ============ ############
UE.samp <- samp[which(samp$SampleGroup=="UE"),]
UE.samprot <- samprot[which(samp$SampleGroup=="UE"),]

PAE.samp <- samp[which(samp$SampleGroup=="PAE"),]
PAE.samprot <- samprot[which(samp$SampleGroup=="PAE"),]

PVE.samp <- samp[which(samp$SampleGroup=="PVE"),]
PVE.samprot <- samprot[which(samp$SampleGroup=="PVE"),]








######### ========= Check that participants have the same order so they can be compared and normalized pairwise ======== #######
(UE.samp$ID == PAE.samp$ID) %>% table()
(VUE.samp$ID == AUE.samp$ID) %>% table()
(PVE.samp$ID == PAE.samp$ID) %>% table() # ohno -> # PVE and PAE includes different patients!!
# Find overlap patient ids between PVE and PAE
PVEPAE.patients <- intersect(PVE.samp$ID, PAE.samp$ID)
(PVE.samp$ID[which(PVE.samp$ID  %in% PVEPAE.patients)] == PAE.samp$ID[which(PAE.samp$ID  %in% PVEPAE.patients)]) %>% table()








############ =============== Calculate differences between pairs of vessels ============= ###########
m.c.diffs <- UE.samprot-PAE.samprot
adj.m.c.diffs <- m.c.diffs
for (patient in 1:nrow(m.c.diffs)){ 
  adj.m.c.diffs[patient,] <- m.c.diffs[patient,] - median(m.c.diffs[patient,]) 
}


marm.c.diffs <- PVE.samprot[which(PVE.samp$ID  %in% PVEPAE.patients),] - PAE.samprot[which(PAE.samp$ID  %in% PVEPAE.patients),]
adj.marm.c.diffs <- marm.c.diffs
for (patient in 1:nrow(marm.c.diffs)){ 
  adj.marm.c.diffs[patient,] <- marm.c.diffs[patient,] - median(marm.c.diffs[patient,]) 
}








############## ============== Write r data file of datasets ============ ############
save(samprot,
     samp,
     prot,
     adj.m.c.diffs,
     adj.marm.c.diffs,
     file="data/lgRFUmedians_placentalclock.RData")
# save 4-vessel data set median normalized venoarterial differences 
