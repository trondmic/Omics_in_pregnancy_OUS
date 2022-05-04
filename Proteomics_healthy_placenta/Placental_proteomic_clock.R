###############
#
# Author: Maren-Helene Langeland Degnes
# Purpose: Use variable selector to find 
#          proteins that can predict gestational age
#
###############


# Load packages
library(tidyverse)
library(pheatmap)
library(glmnet)
library(ggplot2)
library(ggrepel)
require(gridExtra)
library(ROCR)
library(Biobase)
library(c060)
library(parallel)





#####################        Load dataset          ##################
# load median normalized 4-vessel data set, and expression set of STORK data set
load("data/lgRFUmedians.RData")  
load("data/Oslo_STORK_log2win_exprset.RData")
samprot.CSTORK <- exprs(Oslo_STORK_log2win_exprset) %>% t()
samp.CSTORK <- pData(Oslo_STORK_log2win_exprset)
prot <- fData(Oslo_STORK_log2win_exprset)
# Load test restults ~ placental proteome 
load("data/MTT_FBmedians.RData")







###################### Sample validation cohort and run ENSS  ##################
set.seed(123)
val.idx <- which(samp.CSTORK$ID %in% sample(samp.CSTORK$ID %>% unique(),70/3))
samp.tr <- samp.CSTORK[-val.idx,]
tr_y_GA <- samp.tr$GAWeeks
tr_x_mx <- data.frame(samprot.CSTORK[-val.idx,],"BMI"=samp.tr$BMI,"Age"=samp.tr$Age,
                      "Nulliparity"=samp.tr$Nulliparity) %>% as.matrix()
ids <- unique(samp.tr$ID)

samp.val <- samp.CSTORK[val.idx,]
val_y_GA <- samp.val$GAWeeks
val_x_mx <-  data.frame(samprot.CSTORK[val.idx,],"BMI"=samp.val$BMI,"Age"=samp.val$Age,
                        "Nulliparity"=samp.val$Nulliparity)

alphalist <- list()
alphas <- c(1,0.8,0.6,0.4,0.2)

n_trainpart <- nrow(tr_x_mx)
n_valpart <- nrow(val_x_mx)

## Create results container
reslist_within <- list(proteins=list(),
                       models=list(),
                       validate.pred=matrix(ncol=n_valpart,nrow=length(alphas)),
                       validate.truth=matrix(ncol=n_valpart,nrow=length(alphas)))

##### Run ENSS with 5 different alpha values 
for (alpha_i in c(1:5)){
  
  spath <- stabpath(y=tr_y_GA, 
                    x=tr_x_mx, 
                    family="gaussian", 
                    alpha=alphas[alpha_i])
  stabrun <- c060::stabsel(x=spath, error=1, type="pfer", pi_thr=0.6)
  stabproteinsidx <- stabrun$stable
  
  ######### save datasets for linear model which demands me to create these
  tr_dta <- data.frame(y_GA=tr_y_GA,x_mx=tr_x_mx[,stabproteinsidx])
  val_dta <-  data.frame(x_mx=val_x_mx[,stabproteinsidx])
  linearmod <- lm(y_GA~.,data=tr_dta)
  
  ## Save results from elastic net modelling
  reslist_within$proteins[[alpha_i]] <- stabproteinsidx
  reslist_within$models[[alpha_i]] <- linearmod
  reslist_within$validate.pred[alpha_i,] <- predict(linearmod,newdata=val_dta,type="response")
  reslist_within$validate.truth[alpha_i,] <- val_y_GA
  
  ######### perform analysis using LOOCV
  reslist <- list(preds=matrix(ncol=3,nrow=length(ids)),
                  truths=matrix(ncol=3,nrow=length(ids)),
                  proteins=list())
  for (loo_f in 1:length(ids)){
    loo_row <- which(samp.tr$ID==ids[loo_f])
    train_y <- tr_y_GA[-loo_row]
    train_x <- tr_x_mx[-loo_row,]
    
    spath <- stabpath(y=train_y, 
                      x=train_x, 
                      family="gaussian", 
                      alpha=alphas[alpha_i])
    stabrun <- c060::stabsel(x=spath, error=1, type="pfer", pi_thr=0.6)
    stabproteinsidx <- stabrun$stable
    traindta <- data.frame(y_GA=train_y,x_mx=train_x[,stabproteinsidx])
    linearmod <- lm(y_GA~.,data=traindta)
    reslist$proteins[[loo_f]] <- stabproteinsidx
    reslist$preds[loo_f,] <- predict(linearmod,newdata=data.frame(x_mx=tr_x_mx[loo_row,stabproteinsidx]),type="response")
    reslist$truths[loo_f,] <- tr_y_GA[loo_row]
  }
  
  alphalist[[alpha_i]] <- reslist
  
}

save(alphalist,reslist_within,file=paste0("data/06_stability_selection_res_rel_cl_punish_winz.RData"))

















####################### Plot performances of selected alpha ##################

pdf(paste0("results/Figure 5.pdf"),
    width=6.6,height=3.6)
par(mfcol=c(1,2),
    cex.main=0.5,
    cex.axis=0.5,
    cex.lab=0.5,
    mar=c(3,3,2,0.5),
    mgp=c(3, 0.5, 0))
# for (bag_i in 1:3){
load(paste("data/06_stability_selection_res_rel_cl_mandatory_winz.RData"))
for (ai in 1){
  preds <- alphalist[[ai]]$preds %>% t() %>% c()
  truths <- alphalist[[ai]]$truths %>% t() %>% c()
  
  plot(preds,truths,main="LOO CV on training cohort",
       xlab="Predicted GA",ylab="True GA",cex=0.5,lwd=0.5)
  title(xlab="Predicted gestational weeks",line=1.2)
  title(ylab="True gestational weeks",line=1.2)
  abline(a=0,b=1,col="blue")
  lmod <- lm(truths~preds)
  abline(a=lmod$coefficients[1],b=lmod$coefficients[2],col="orange")
  legend("bottomright",legend=paste("Pearson correlation=",cor(preds,truths) %>% round(.,3)),bty="n",cex=0.4)
  legend("topleft",pch=c(1,NA,NA),lty=c(NA,1,1),
         legend=c("GA in weeks","Curve from a perfect prediction","Actual curve estimated based on lm"),
         col=c(1,"blue","orange"),bty="n",cex=0.4)
  
  text(7,36,"A.",xpd=T,cex=0.5)
  
  plot(reslist_within$validate.pred[ai,],
       reslist_within$validate.truth[ai,],
       xlab="Predicted GA",ylab="True GA",cex=0.5,lwd=0.5,
       main="Trained model on validation cohort",
       xlim=c(10,36),
       ylim=c(12,34))
  title(xlab="Predicted gestational weeks",line=1.2)
  title(ylab="True gestational weeks",line=1.2)
  abline(a=0,b=1,col="blue")
  val.preds <- reslist_within$validate.pred[ai,]
  val.truths <- reslist_within$validate.truth[ai,]
  lmod <- lm(val.truths~val.preds)
  abline(a=lmod$coefficients[1],b=lmod$coefficients[2],col="orange")
  legend("bottomright",legend=paste("Pearson correlation=",cor(reslist_within$validate.pred[ai,],
                                                               reslist_within$validate.truth[ai,]) %>% round(.,3)),bty="n",cex=0.4)
  legend("topleft",pch=c(1,NA,NA),lty=c(NA,1,1),
         legend=c("GA in weeks","Curve from a perfect prediction","Actual curve estimated based on lm"),
         col=c(1,"blue","orange"),bty="n",cex=0.5)
  text(7,36,"B.",xpd=T,cex=0.5)
  
};dev.off()













