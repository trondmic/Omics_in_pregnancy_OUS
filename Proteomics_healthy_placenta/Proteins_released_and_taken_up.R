##############
#
#
# Author: Maren-Helene Langeland Degnes
# Purpose: Test venoarterial differences and 
#           define proteins released and taken up
#           by the placenta
#
##############







# Load packages
library(haven)
library(tidyverse)
library(scales)
library(openxlsx)





options(stringsAsFactors = FALSE)
rm(list = ls())





###################### ================ Function for Welch's t-tests ============ ###############
PairMulTestD <- function(difftab,difftab2=NA){
  test <- c()
  if(length(difftab2)<2){
    for (i in 1:ncol(difftab)){
      each.test <- t.test(difftab[,i])
      test <- rbind(test, c(each.test$estimate, 
                            each.test$statistic, 
                            each.test$conf.int[1], 
                            each.test$conf.int[2], 
                            each.test$p.value,NA))
    }
    test[,6] <- p.adjust(test[,5], method="BH")
    colnames(test) <- c("mean of differences", 
                        "t statistic", 
                        "lower .95 CI", 
                        "upper .95 CI", 
                        "p value", 
                        "adjusted p")
  } else{  
    
    for (i in 1:ncol(difftab)){
    each.test <- t.test(difftab[,i],difftab2[,i])
    test <- rbind(test, c(each.test$estimate[1],
                          each.test$estimate[2],
                          each.test$statistic, 
                          each.test$conf.int[1], 
                          each.test$conf.int[2], 
                          each.test$p.value,NA))
    }
    test[,7] <- p.adjust(test[,6], method="BH")
    colnames(test) <- c("mean of differences1", 
                        "mean of differences1", 
                        "t statistic", 
                        "lower .95 CI", 
                        "upper .95 CI", 
                        "p value", 
                        "adjusted p")}
  test <- test %>% as.data.frame()
  
  
  return(test)
}

  





############## =============== Perform and save results from tests ================ ##############
# load 4-vessel data set median normalized venoarterial differences 
load("data/lgRFUmedians_placentalclock.RData")

MatNormTest <- PairMulTestD(difftab = adj.m.c.diffs)
ArmNormTest <- PairMulTestD(difftab = adj.marm.c.diffs)

secrprot <- which(MatNormTest$`adjusted p`<0.05 & MatNormTest$`t statistic`>0)

MatPLARMtest <- PairMulTestD(difftab = adj.m.c.diffs[,secrprot], 
                             difftab2 = adj.marm.c.diffs[,secrprot])
MatNormTest$PLARMpadj <- NA
MatNormTest$PLARMt <- NA
MatNormTest$PLARMpadj[secrprot] <- MatPLARMtest$`adjusted p`
MatNormTest$PLARMt[secrprot] <- MatPLARMtest$`t statistic`

save(MatNormTest,
     ArmNormTest,
     MatPLARMtest,
     file="data/MTT_FBmedians_placentanclock.RData")
# Save test results as MTT








############ ====================== p value histogram ============= ############

pdf('results/Supplementary figure 2.pdf',
     width=3.3, height=3)
load(paste0("data/MTT_FBmedians_placentanclock.RData"))
b <- 0.025
ylims <- c(0,1000)
par(mfcol=c(1,1),
    mar=c(2,2,1,0.5),
    cex.axis=0.5,
    mgp=c(3,0.5,0),
    cex.lab=0.5,
    cex.main=0.8)

p <- rep(NA, length(MatNormTest$`p value`))
hist(MatNormTest$`p value`, 
     breaks=seq(0,1,b), 
     main=paste0(""),
     xlab="P-values",
     ylim = ylims)
title(ylab="Frequency",line=1.2)
title(xlab="P-values",line=1)
dev.off()




