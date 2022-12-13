######
#
# Author: Maren-Helene Langeland Degnes
# Date: 24.Oct.22
#
####
library(tidyverse)
library(splines)
library(lme4)
library(Biobase)

# different values of knots placements to loop over
innerknots <-c(NA,24.5,12.1,14.1,15.1)
firstknots <- c(12,12,12,14,15)
lastknots <- c(34,34,34,34,34)
prins <- c(1419)
mains <- c("A. Only boundary knots",
           "B. Inner knot in middle of boundary knots",
           "C. Inner knot equal to first boundary knot",
           "D. Inner knot .1 more than first boundary knot",
           "E. First boundary + inner knot placed at later GW")



par(mfcol=c(length(prins)+1,length(innerknots)),mar=c(4,4,4,1),cex.lab=1,cex.main=1)
load("STORKdataset.RData")
dtamx <- exprs(STORKdataset) %>% t()
ano <- pData(STORKdataset)
prot <- fData(STORKdataset)
sub.idx <- which(ano$PEgroups %in% c(0,2))  # only late-onset preeclampsia
sub.ano <- ano[sub.idx,]
sub.dtamx <- dtamx[sub.idx,]


GA <- sub.ano$GAWeeks %>% sort()

# Plot protein trajectories
for (k in 1:length(firstknots)){
  
  #alternativerly replace naturalsplinesmx with bsplinesmx
  bsplinesmx <- bs(GA,
                   knot=innerknots[k],
                   degree=2,
                   Boundary.knots=c(firstknots[k],lastknots[k]))
  plot(c(),xlim=range(GA),ylim=c(-1.5,1.5),
       main=mains[k],
       ylab="f(gestational weeks)",xlab="Gestational weeks")
  matlines(GA, bsplinesmx)
  lines(GA, rowSums(bsplinesmx), col = "purple",lwd=2)
  abline(v=c(firstknots[k],innerknots[k],lastknots[k]),col=4,lty=2)
  text(firstknots[k]+1,1.5, firstknots[k],col=4)
  fi.knot.diff <-ifelse(k%in%c(1:2),1,2)
  text(innerknots[k]+fi.knot.diff,1.5, innerknots[k],col=4,adj=0)
  text(lastknots[k]-1,1.5, lastknots[k],col=4)
    
  for ( pr in prins){
  
  library(tidyverse)
  library(lme4)
  library(splines)
  library(Biobase)
  source("functions.R")  # My function patientplots() to plot the protein abundance per participant across gestational weeks

    
    ######################### ===================== Load data, fit mixed model to each protein and test estimated effects ======= #################
    ########### PRE ANALYSIS
    ## Load data
    load("STORKdataset.RData")
    dtamx <- exprs(STORKdataset) %>% t()
    ano <- pData(STORKdataset)
    prot <- fData(STORKdataset)
    sub.idx <- which(ano$PEgroups %in% c(0,2))  # only late-onset preeclampsia
    sub.ano <- ano[sub.idx,]
    sub.dtamx <- dtamx[sub.idx,]
    
    
    # Create containers needed inside for-loop
    ST.pred.mx <- c()
    pred.mx <- c()
    models <- list()
    modelsi <- 1
    bsreslist <- list()   # container for testing coefficients in mm
        
    # Save objects for mixed model
    ID <- sub.ano$ID
    PE <- sub.ano$PE %>% as.numeric() %>% as.factor()
    BMI <- sub.ano$BMI %>% scale()           #standardize BMI
    age <- sub.ano$Age %>% as.numeric() %>% scale()  #standardize age
    nulliparity <- sub.ano$Nulliparity %>% as.factor()
    meanTP <- tapply(sub.ano$GAWeeks %>% as.numeric(),sub.ano$TimePoint,mean) 
    prot.indeces <- 1:ncol(sub.dtamx)
    
    
    # par(mfrow=c(1,1))
    # pr <- 4008
    y <- sub.dtamx[,pr]
    x <- sub.ano$GAWeeks
    
    
    ########## test all kinds of settings to see which gives hits ##########
    knot <- innerknots[k]
    b.knot <- c(firstknots[k],lastknots[k])
    xnew1 <- seq(firstknots[k],lastknots[k],length.out=100)
    bsplines <- bs(x, degree=2,knots=c(knot), Boundary.knots =c(firstknots[k],lastknots[k]))
    bsplines.df <- bsplines %>% as.data.frame()
    colnames(bsplines.df) <- paste0("bs",1:ncol(bsplines.df))
    design <- cbind(bsplines.df,ID,PE,y,BMI,age)
    
    if(ncol(bsplines.df)==2){
      mmfit_var1 <- lmer(y ~ 1 + bs1*PE + bs2*PE + BMI + age + (1| ID), data=design,REML=F)
    }
    if(ncol(bsplines.df)==3){
      mmfit_var1 <- lmer(y ~ 1 + bs1*PE + bs2*PE + bs3*PE+ BMI + age + (1| ID), data=design,REML=F)
    } 
    
    lmerTest::as_lmerModLmerTest(mmfit_var1) %>% summary()
    
    # Predict curves
    predsplines_var1 <- predict(bsplines,xnew1)
    colnames(predsplines_var1) <- paste0("bs",1:ncol(bsplines.df))
    
    patientplots(pr,PE=T,earlyPE4V=F,main=prot$EntrezGeneSymbol[pr])
    # legend("topleft",legend=paste0("Knots=",firstknots[k],",",innerknots[k],",",lastknots[k]),bty="n")
    
    BMI <- rep(mean(BMI),100)
    age <- rep(mean(age),100)
    
    PE <- rep(0,length(xnew1)) %>% as.factor()
    predcurve_C_var1 <- predict(mmfit_var1,predsplines_var1, re.form=~0)
    PE <- rep(1,length(xnew1)) %>% as.factor()
    predcurve_PE_var1 <- predict(mmfit_var1,predsplines_var1, re.form=~0)
    lines(xnew1,predcurve_C_var1,lwd=2,col="black")
    lines(xnew1,predcurve_PE_var1,col="red",lwd=2)  
  } 
}

