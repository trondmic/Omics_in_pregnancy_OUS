#########
#
# Author: Ina Jungersen Andresen
# Date: 2.dec.2022
#
# Purpose: Test if proteins are differential abundant per visit using LIMMA tests
#
#
######


# load packages
library(tidyverse)
library(writexl)
library(readxl)
library(Biobase)
library(limma)
library(ggVennDiagram)
library(pcaMethods)
library(scales)





#Make model for Limma test, normal vs preeclampsia at different time points. Adjust for age, nulliparity and BMI --------
load("data/MOM_STORK_log2win_exprset.RData")

#Remove PE-early samples
MOM_STORK_log2win_exprset <- MOM_STORK_log2win_exprset[, pData(MOM_STORK_log2win_exprset)$PEgroups != "1"]

#Change the 1, 2, 3 to T1, T2, T3 in the time column
MOM_STORK_log2win_exprset$TimePoint[MOM_STORK_log2win_exprset$TimePoint == "1"] <- "T1"
MOM_STORK_log2win_exprset$TimePoint[MOM_STORK_log2win_exprset$TimePoint == "2"] <- "T2"
MOM_STORK_log2win_exprset$TimePoint[MOM_STORK_log2win_exprset$TimePoint == "3"] <- "T3"

#Merge health (preeclamptic or normal) with Time
Time_health <- as.factor(paste(MOM_STORK_log2win_exprset$PE, MOM_STORK_log2win_exprset$TimePoint, sep = "."))

#Extract the sample information
samp_STORK <- pData(MOM_STORK_log2win_exprset)

#Make the design matrix and adjust for age, parity, smoking and BMI
design <- model.matrix(~0+Time_health+as.numeric(Age)+as.factor(Nulliparity)+as.numeric(BMI),samp_STORK)
colnames(design)[colnames(design) %in% c("Time_health0.T1", "Time_health0.T2",
                                         "Time_health0.T3", "Time_health1.T1",
                                         "Time_health1.T2", "Time_health1.T3", 
                                         "as.numeric(Age)", "as.factor(Nulliparity)1",
                                         "as.numeric(BMI)")] <- c("Control.T1", "Control.T2", 
                                                                  "Control.T3", "PE.T1", 
                                                                  "PE.T2", "PE.T3", 
                                                                  "Age", "Nulliparity", 
                                                                  "BMI")


#See this about blocking: https://mailman.stat.ethz.ch/pipermail/bioconductor/2014-February/057887.html
#Estimate the correlation between time measurements made on the same patient
#Treats patients (ID) as random effect
corfit <- duplicateCorrelation(MOM_STORK_log2win_exprset, design, block = MOM_STORK_log2win_exprset$ID)
corfit$consensus 

# Limma t-test ------------------------------------------------------------
#Do the limma test, and block for correlations effects from patients
fit <- lmFit(MOM_STORK_log2win_exprset, design, block = MOM_STORK_log2win_exprset$ID, correlation = corfit$consensus)

#Make the contrast: Which proteins are deferentially expressed between preeclamptic and and normal patients 
#at the different time point? 
contrast <- makeContrasts(T1 = PE.T1 - Control.T1, 
                          T2 = PE.T2 - Control.T2, 
                          T3 = PE.T3 - Control.T3,
                          levels = design)


#Compute the moderated t-tests
fit <- contrasts.fit(fit, contrast)
efit <- eBayes(fit)

save(fit, design, efit, file = "fit_STORK_MoM_adjustedAgeParityBMI_blockedIDs_PElate.RData")


#Extract results -------------------------------------------------
load("fit_STORK_MoM_adjustedAgeParityBMI_blockedIDs_PElate.RData")

#Negative fold change means higher in controls
#Positive fold change means higher in PE

#No cut off values
T1 <- topTable(efit, coef = "T1", number = "all", adjust.method="BH")
T2 <- topTable(efit, coef = "T2", number = "all", adjust.method="BH")
T3 <- topTable(efit, coef = "T3", number = "all", adjust.method="BH")

T1 <- T1 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


T2 <- T2 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


T3 <- T3 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T1, path = "T1.xlsx",
           col_names = TRUE)
write_xlsx(T2, path = "T2.xlsx",
           col_names = TRUE)
write_xlsx(T3, path = "T3.xlsx",
           col_names = TRUE)


#Cutoff FDR < 0.05
#Cut-off values: q-value < 0.05 
results <- decideTests(efit, adjust.method = "BH", p.value = 0.05)
vennDiagram(results, main = "FDR<0.05")

T1 <- topTable(efit, coef = "T1", number = "all", adjust.method="BH", p.value = 0.05)
T2 <- topTable(efit, coef = "T2", number = "all", adjust.method="BH", p.value = 0.05)
T3 <- topTable(efit, coef = "T3", number = "all", adjust.method="BH", p.value = 0.05)

T1 <- T1  %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T2 <- T2  %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T3 <- T3  %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T1, path = "T1_FDR0.05.xlsx",
           col_names = TRUE)
write_xlsx(T2, path = "T2_FDR0.05.xlsx",
           col_names = TRUE)
write_xlsx(T3, path = "T3_FDR0.05.xlsx",
           col_names = TRUE)


#Unique proteins:
T3_proteins <- T3$EntrezGeneSymbol
T3_proteins_unique <- unique(T3_proteins)

T3_proteins_higher_PE <- T3 %>%
  filter(logFC > 0)
T3_proteins_higher_PE_unique <- unique(T3_proteins_higher_PE$EntrezGeneSymbol)

T3_proteins_lower_PE <- T3 %>%
  filter(logFC < 0)
T3_proteins_lower_PE_unique <- unique(T3_proteins_lower_PE$EntrezGeneSymbol)


#Cut off values: q-value = 0.25 and log2FC = 0.5 
results <- decideTests(efit, adjust.method = "BH", p.value = 0.25, lfc = 0.5)
vennDiagram(results, main = "FDR<0.25, logFC >0.5")

T1 <- topTable(efit, coef = "T1", number = "all", adjust.method="BH", p.value = 0.25, lfc = 0.5)
T2 <- topTable(efit, coef = "T2", number = "all", adjust.method="BH", p.value = 0.25, lfc = 0.5)
T3 <- topTable(efit, coef = "T3", number = "all", adjust.method="BH", p.value = 0.25, lfc = 0.5)

T1 <- T1 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T2 <- T2 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T3 <- T3  %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

write_xlsx(T1, path = "T1_FDR0.25_log2fc0.5.xlsx",
           col_names = TRUE)
write_xlsx(T2, path = "T2_FDR0.25_log2fc0.5.xlsx",
           col_names = TRUE)
write_xlsx(T3, path = "T3_FDR0.25_log2fc0.5.xlsx",
           col_names = TRUE)



