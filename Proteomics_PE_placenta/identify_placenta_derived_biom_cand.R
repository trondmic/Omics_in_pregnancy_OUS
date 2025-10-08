#################

#
# Author: Ina Jungersen Andresen
# Date: 11.04.25
# Purpose: Using moderated t-tests (limma R-package version 3.50.0) 37  , 
#          we compared the antecubital vein log2 RFUs of the placenta-derived 
#          proteins at each of the three visits, well before the onset of PE, 
#          between PE and controls in the STORK cohort adjusting for the known 
#          risk factors; maternal age, BMI and nulliparity. 
#
#
#
#################
 

library(tidyverse)
library(writexl)
library(readxl)
library(Biobase)
library(limma)
library(ggVennDiagram)
library(scales)
library(statmod)
library(gridExtra)
library(ggrepel)
library(grid)
 
 
# Adjust for Age, nulliparity, and BMI ------------------------------------
load("data/exprset.RData")
#load list of placenta-derived proteins
released_proteins <- read_xlsx("Significantly_released_proteins.xlsx",
                               sheet = 1, col_names = FALSE)
released_proteins <- released_proteins$...1
 
#Filter the expressionset on the 620 placenta-derived proteins 
lgRFU_released_proteins <- exprset[rownames(exprset) %in% released_proteins,]
 
#Continue with only STORK samples 
lgRFU_released_proteins <- lgRFU_released_proteins[, (pData(lgRFU_released_proteins)$STORK %in% "TRUE")]
 
#Remove the patients with early PE
lgRFU_released_proteins <- lgRFU_released_proteins[, !(pData(lgRFU_released_proteins)$PEgroups %in% "1")]
samp <- pData(lgRFU_released_proteins)
lgRFU <- exprs(lgRFU_released_proteins)
 
#Check numbers
table(samp$TimePoint)
table(samp$PEgroups)
 
#Change the 1, 2, 3 to T1, T2, T3 in the time column
lgRFU_released_proteins$TimePoint[lgRFU_released_proteins$TimePoint == "1"] <- "T1"
lgRFU_released_proteins$TimePoint[lgRFU_released_proteins$TimePoint == "2"] <- "T2"
lgRFU_released_proteins$TimePoint[lgRFU_released_proteins$TimePoint == "3"] <- "T3"
 
#Merge health (preeclamptic or normal) with Time
Time_health <- as.factor(paste(lgRFU_released_proteins$PE, lgRFU_released_proteins$TimePoint, sep = "."))
 
#Make the design matrix and adjust for age, parity, smoking and BMI
design <- model.matrix(~0+Time_health+as.numeric(Age)+as.factor(Nulliparity)+as.numeric(BMI),samp)
colnames(design)[colnames(design) %in% c("Time_health0.T1", "Time_health0.T2",
                                         "Time_health0.T3", "Time_health1.T1",
                                         "Time_health1.T2", "Time_health1.T3", 
                                         "as.numeric(Age)", "as.factor(Nulliparity)1",
                                         "as.numeric(BMI)")] <- c("Control.T1", "Control.T2", 
                                                                  "Control.T3", "PE.T1", 
                                                                  "PE.T2", "PE.T3", 
                                                                  "Age", "Nulliparity", 
                                                                  "BMI")
 
 
#Estimate the correlation between time measurements made on the same patient
#Treats patients (ID) as random effect
corfit <- duplicateCorrelation(lgRFU_released_proteins, design, block = lgRFU_released_proteins$ID)
corfit$consensus 
#0.4266294 Adjusted for Age, Nulliparity and BMI, PE-late
 
#Do the limma test, and block for correlations effects from patients
fit <- lmFit(lgRFU_released_proteins, design, block = lgRFU_released_proteins$ID, correlation = corfit$consensus)
 
#Make the contrast: Which proteins are deferentially expressed between preeclamptic and and normal patients 
#at the different time point? 
contrast <- makeContrasts(T1 = PE.T1 - Control.T1, 
                          T2 = PE.T2 - Control.T2, 
                          T3 = PE.T3 - Control.T3,
                          levels = design)
 
 
#Compute the moderated t-tests
fit <- contrasts.fit(fit, contrast)
efit <- eBayes(fit)
 
 
#Negative fold change means higher in PE
#Positive fold change means higher in Controls
 
#Print results and filter on adjusted p-value 0.05
T1 <- topTable(efit, coef = "T1", number = "all", p.value = 0.05, adjust.method="BH")
T2 <- topTable(efit, coef = "T2", number = "all", p.value = 0.05, adjust.method="BH")
T3 <- topTable(efit, coef = "T3", number = "all", p.value = 0.05, adjust.method="BH")
 
T1 <- T1 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
 
 
T2 <- T2 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
 
 
T3 <- T3 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
