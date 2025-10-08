####### script information ##########
#
# Author: Maren-Helene Langeland Degnes
# Date: 17.12.2021
# Purpose: A. LIMMA tests across 5000 proteins for controls and PE to identify protein released or taken up by the placenta
#
#
############



#### load packages ########
options(stringsAsFactors = FALSE)
library(haven)
library(tidyverse)
library(scales)
library(openxlsx)
library(limma)
library(readxl)
library(gridExtra)
library(ggrepel)
library(writexl)
library(grid)


######  Load data made in preprocessing.R ######
pth <- "" # path to data
load(paste0(pth,"data/lgRFUmedianswins.RData"))







####### A. Identify proteins taken up and/or released by the placenta in each group #####

#### Test pl differences in C samples
UE.samp <- samp[c(which(samp$PE==0 & samp$SampleGroup=="UE")),]
table(rownames(UE.samp)== c(rownames(adj.m.c.diffs)))
design <- model.matrix(~1,UE.samp)
colnames(design) <- c("Intercept")
C_matp_limma_cl <- adj.m.c.diffs %>% 
  t() %>% 
  lmFit(.,design=design) %>% 
  eBayes() %>% 
  topTable(.,coef="Intercept",number="all", adjust.method = "BH",sort.by="none")

#### Test pl differences in PE samples
UE.samp <- samp[c(which(samp$PE==1 & samp$SampleGroup=="UE")),]
table(rownames(UE.samp)== c(rownames(adj.m.pe.diffs)))
design <- model.matrix(~1,UE.samp)
colnames(design) <- c("Intercept")
PE_matp_limma_cl <- adj.m.pe.diffs %>% 
  t() %>% 
  lmFit(.,design=design) %>% 
  eBayes() %>% 
  topTable(.,coef="Intercept",number="all", adjust.method = "BH",sort.by="none")


restab <-   data.frame(Proteinname = prot$TargetFullName, 
                       GeneSymbol = prot$EntrezGeneSymbol,
                       SomaID = prot$SomaId,
                       AptamerName = prot$AptName,
                       "Log2FC_C" = C_matp_limma_cl$logFC,
                       "RFUrati_C" = apply(2^adj.m.c.diffs,2,mean),
                       "BHAdjP_C" = C_matp_limma_cl$adj.P.Val,
                       "Log2FC_PE" = PE_matp_limma_cl$logFC,
                       "RFUrati_PE" = apply(2^adj.m.pe.diffs,2,mean),
                       "BHAdjP_PE" = PE_matp_limma_cl$adj.P.Val)
write.xlsx(restab,file=paste0(pth,"restable.xlsx"))




