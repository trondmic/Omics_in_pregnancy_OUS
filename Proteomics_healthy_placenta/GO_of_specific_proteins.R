###############
#
#
# Author: Ina Jungersen Andresen
# Purpose: Run gene ontology enrichment tests on 
#           placenta specific proteins
#
#############

library(org.Hs.eg.db)
library(readxl)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)

#Create background list containing EntrezGeneIDs and protein name (here AptName) for all the measured proteins
load("data/MOM_Oslo_STORK_log2win_exprset.RData")
background <- fData(MOM_Oslo_STORK_log2win_exprset)
background <- background %>%
  dplyr::select(TargetFullName, AptName, EntrezGeneID, EntrezGeneSymbol, SomaId, UniProt, Target)
rm(MOM_Oslo_STORK_log2win_exprset)

background_entrez <- background %>%
  select(EntrezGeneID, AptName)

#Upload the list of proteins you want to de the enrichment analysis on
upreg_spec <- read_xlsx("data/enrichment.xlsx", sheet = 3, col_names = FALSE)
upreg_spec <- upreg_spec[,c(1:4)]
colnames(upreg_spec) <- c("TargetFullName", "UniProt", "SomaId", "AptName")

#Add EntrezGeneIds to the list of proteins by merging with backround list
upreg_spec <- merge(upreg_spec, background_entrez, by = "AptName")

#Do an over representation analysis (ORA) with the enrichGO function from clusterProfiler
#Extract GO-terms belong to "Biological process" only
#Use the human genome as reference
upreg_spec_enrichgo_BP <- enrichGO(upreg_spec$EntrezGeneID, OrgDb = "org.Hs.eg.db", ont = "BP", universe = background$EntrezGeneID, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)

#Simplify the GO-terms
upreg_spec_enrichgo_BP_simp <- simplify(upreg_spec_enrichgo_BP) 

#Calculate rich factor 
upreg_spec_enrichgo_BP_simp <- mutate(upreg_spec_enrichgo_BP_simp, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) 

#Filter the data to contain enriched GO terms with rich factor > 0.15
BP_simp_filt <- upreg_spec_enrichgo_BP_simp %>%
  filter(richFactor > 0.15)

#Plot a dotplot
dotplot_BP_simp <- dotplot(BP_simp_filt, x = "richFactor", showCategory = 100) 

#Save the plot
ggsave(filename= "results/Figure 2.pdf",
       plot = dotplot_BP_simp,
       device = "pdf",
       path = "data/",
       width = 20,
       height = 30,
       units = "cm")


