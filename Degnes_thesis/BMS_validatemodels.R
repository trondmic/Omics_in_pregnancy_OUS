
#######
#
# Author: Maren-Helene Langeland Degnes
# Date: 13. Dec. 2022
# Purpose: Obtain model performances of the statistical machine learning models 
#          elastic net and random forest, and the elastic net with stability selection
#
#########

require(tidyverse)
require(randomForest)
require(ROCR)
require(glmnet)
require(Biobase)
require(parallel)
require(c060)
source("~/c060_functions.R")



############### LOAD data and sample train and test###############
load("detrendedSTORKdataset.RData")
dta <- exprs(STORKdataset) %>% t()  # proteins in columns and samples in rows
ano <- pData(STORKdataset)
prot <- fData(STORKdataset)
colnames(dta) <- prot$AptName

dta <- dta[which(ano$PEgroups %in% c(0,2)),]
ano <- ano[which(ano$PEgroups %in% c(0,2)),]
ranges <- tapply(ano$GAWeeks , ano$TimePoint, range)
table(rownames(MoMmx) == rownames(ano))        # Check

###  sample FOLDIDs
ano$FOLDID <- NA
all_ids <- ano$ID %>% unique()
set.seed(1,kind="default")                            # set seed so the sampling is reproducible 
for ( id in all_ids ){  ano$FOLDID[which(ano$ID == id)] <- sample(1:5,1)}
table(ano$FOLDID,ano$PEgroups)


nseeds <- 1





######################### Fit all models and collect performances #########################
for ( tp in 1:3){
  
  # Choose data on time point tp and only predrols and late PE
  tp.latePE.idx <- which(ano$TimePoint==tp)
  y_mx <- ano[tp.latePE.idx,]
  x_mx <- dta[tp.latePE.idx,]
  
  # Check
  table(rownames(x_mx) == rownames(y_mx))
  
  # Set clinical data and number of folds for CV
  x <- data.frame(x_mx,BMI=y_mx$BMI,age=as.numeric(y_mx$Age),parity=y_mx$Nulliparity%>%as.factor())
  x_nocl <- data.frame(x_mx)
  y <- y_mx$PE %>% as.factor()
  n <- length(y)
  K <- 5
  
  for ( alph in c(1,0.8,0.6,0.4,0.2) ){
    
    rf.list <- list(); for (i in nseeds){rf.list[[i]] <- list()}
    EN.list <- list(); for (i in c(nseeds,nseeds+length(nseeds))){EN.list[[i]] <- list()}
    ENSS.list <- list(); for (i in c(nseeds,nseeds+length(nseeds))){ENSS.list[[i]] <- list()}
    ENSS.list.ratio <- list(); for (i in c(nseeds,nseeds+length(nseeds))){ENSS.list.ratio[[i]] <- list()}
    ENSS.list.ratiomand <- list(); for (i in c(nseeds,nseeds+length(nseeds))){ENSS.list.ratiomand[[i]] <- list()}
    logcl.list.ratio <- list(); for (i in nseeds){logcl.list.ratio[[i]] <- list()}
    logcl.list.bms <- list(); for (i in nseeds){logcl.list.bms[[i]] <- list()}
    logcl.list.ratio.bms <- list(); for (i in nseeds){logcl.list.ratio.bms[[i]] <- list()}
    K.y <- list(); for (i in nseeds){K.y[[i]] <- list()}
    
    for (nseed in nseeds){
      
      # Run training and testing
      for (k in 1:K){
        
        x.tr <- x[-which(y_mx[,47+nseed] == k),] %>% as.data.frame()
        x.te <- x[which(y_mx[,47+nseed] == k),] %>% as.data.frame()
        
        x.tr.nocl <- x_nocl[-which(y_mx[,47+nseed] == k),] %>% as.data.frame()
        x.te.nocl <- x_nocl[which(y_mx[,47+nseed] == k),] %>% as.data.frame()
        
        y.tr <- y[-which(y_mx[,47+nseed] == k)] %>% as.factor()
        y.te <- y[which(y_mx[,47+nseed] == k)] %>% as.factor()
        K.y[[nseed]][[k]] <- y.te
        
        
        
        if(alph == 1){            #Run random forests only once, not for every alpha
          
          # i Random forests!
          set.seed(500)
          rf.mod <- randomForest(x= x.tr,  y= y.tr,  ntree=1000, keep.forest=TRUE, importance=TRUE)
          rf.list[[nseed]][[k]] <- predict(rf.mod,  newdata= x.te,  type="prob")
          
          
          # ii Logistic regression including FLT1/PGF ratio + clinical variables
          ratio <- x.tr[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.tr[,2505]
          clinical_variables <- x.tr[,which(colnames(x.tr)%in%c("BMI","parity","age"))]
          train_dataset <- cbind(ratio,clinical_variables)
          log.model <- glm(y.tr ~ .,  data=train_dataset, family = binomial)
          ratio <- x.te[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.te[,2505]
          clinical_variables <- x.te[,which(colnames(x.te)%in%c("BMI","parity","age"))]
          test_dataset <- cbind(ratio, clinical_variables)
          logcl.list.ratio[[nseed]][[k]] <- predict(log.model, newdata = test_dataset, type="response")
          
          
          # iii Logistic regression including the biomarker signature
          clinical_variables <- x.tr[,which(colnames(x.tr)%in%c("BMI","parity","age"))]
          biomarker_signature_proteins <- x.tr[,list(c(),c(1319,3159,3797),c(1480, 1544, 3159, 3928, 4008))[[tp]]]
          train_dataset <- cbind(biomarker_signature_proteins,clinical_variables)
          log.model <- glm(y.tr ~ .,  data=train_dataset, family = binomial)
          clinical_variables <- x.te[,which(colnames(x.te)%in%c("BMI","parity","age"))]
          biomarker_signature_proteins <- x.te[,list(c(),c(1319,3159,3797),c(1480, 1544, 3159, 3928, 4008))[[tp]]]
          test_dataset <- cbind(biomarker_signature_proteins, clinical_variables)
          logcl.list.bms[[nseed]][[k]] <- predict(log.model, newdata = test_dataset, type="response")
          
          
          # iv Logistic regression including the biomarker signature + FLT1/PGF ratio + clinical variables
          ratio <- x.tr[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.tr[,2505]
          clinical_variables <- x.tr[,which(colnames(x.tr)%in%c("BMI","parity","age"))]
          biomarker_signature_proteins <- x.tr[,list(c(),c(1319,3159,3797),c(1480, 1544, 3159, 3928, 4008))[[tp]]]
          train_dataset <- cbind(ratio,biomarker_signature_proteins,clinical_variables)
          log.model <- glm(y.tr ~ .,  data=train_dataset, family = binomial)
          ratio <- x.te[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.te[,2505]
          clinical_variables <- x.te[,which(colnames(x.te)%in%c("BMI","parity","age"))]
          biomarker_signature_proteins <- x.te[,list(c(),c(1319,3159,3797),c(1480, 1544, 3159, 3928, 4008))[[tp]]]
          test_dataset <- cbind(ratio, biomarker_signature_proteins, clinical_variables)
          logcl.list.ratio.bms[[nseed]][[k]] <- predict(log.model, newdata = test_dataset, type="response")
          
        }
        
        set.seed(500)
        # v Run elastic net including LASSO
        EN.mod <- cv.glmnet(x=x.tr %>% data.matrix(),y=y.tr,alpha=alph,type.measure="class",
                            family="binomial",penalty.factor=c(rep(1,ncol(MoMmx)),rep(0,3)))
        EN.list[[nseed]][[k]] <- predict(EN.mod$glmnet.fit,s=EN.mod$lambda.min,
                                          newx=x.te%>%data.matrix(), type="response")
        EN.coef <- predict(EN.mod, s = EN.mod$lambda.min , type = "coefficients")
        EN.coef <- as.matrix(EN.coef)                       # Prepare the prediction for matrix-indexing
        EN.list[[nseed+length(nseeds)]][[k]] <- EN.coef[EN.coef != 0,]
        
        
        # vi Run elastic net with stability selection (ENSS) and then logistic regression
        set.seed(500)
        source("K:/Sensitivt/Forskning/2012-5678_Placenta/Proteomics/PredictPreeclampsia/scripts/04_c060_functions.R")
        spath <- stabpath(y=y.tr, x=x.tr, size = 0.632,
                          steps = 100,mc.cores=1, family="binomial",
                          weakness=.8, alpha=alph)
        stabproteinsidx <- c060::stabsel(spath, error=1, type="pfer", pi_thr=0.6)$stable
        names(stabproteinsidx) <- c(prot$AptName,"BMI","Age","Nulliparity")[stabproteinsidx]
        log.model <- glm(y.tr ~ .,  data=x.tr[,stabproteinsidx], family = binomial)
        ENSS.list[[nseed]][[k]] <- predict(log.model, newdata = x.te[,stabproteinsidx], type="response")
        ENSS.list[[nseed+length(nseeds)]][[k]] <- stabproteinsidx
        
        
        # vii Run elastic net with stability selection and then logistic regression incl ratio in input matrix
        set.seed(500)
        source("K:/Sensitivt/Forskning/2012-5678_Placenta/Proteomics/PredictPreeclampsia/scripts/04_c060_functions.R")
        ratio <- x.tr[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.tr[,2505]
        clinical_variables <- x.tr[,which(colnames(x.tr)%in%c("BMI","parity","age"))]
        proteins <- x.tr[,which(!colnames(x.tr)%in%c("BMI","parity","age"))]
        train_dataset <- cbind(proteins,ratio,clinical_variables)
        spath <- stabpath(y=y.tr, x=train_dataset, size = 0.632,
                          steps = 100,mc.cores=1, family="binomial",
                          weakness=.8, alpha=alph)
        stabproteinsidx <- c060::stabsel(spath, error=1, type="pfer", pi_thr=0.6)$stable
        names(stabproteinsidx) <- c(prot$AptName,"BMI","Age","Nulliparity")[stabproteinsidx]
        log.model <- glm(y.tr ~ .,  data=train_dataset[,stabproteinsidx], family = binomial)
        ratio <- x.te[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.te[,2505]
        clinical_variables <- x.te[,which(colnames(x.te)%in%c("BMI","parity","age"))]
        proteins <- x.te[,which(!colnames(x.te)%in%c("BMI","parity","age"))]
        test_dataset <- cbind(proteins,ratio,clinical_variables)
        ENSS.list.ratio[[nseed]][[k]] <- predict(log.model, newdata = test_dataset[,stabproteinsidx], type="response")
        ENSS.list.ratio[[nseed+length(nseeds)]][[k]] <- stabproteinsidx
        
        
        # viii Run elastic net with stability selection and then logistic regression incl ratio as mandatory variable
        set.seed(500)
        source("K:/Sensitivt/Forskning/2012-5678_Placenta/Proteomics/PredictPreeclampsia/scripts/04_c060_functions_4lastmandatory.R")
        ratio <- x.tr[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.tr[,2505]
        clinical_variables <- x.tr[,which(colnames(x.tr)%in%c("BMI","parity","age"))]
        proteins <- x.tr[,which(!colnames(x.tr)%in%c("BMI","parity","age"))]
        train_dataset <- cbind(proteins,ratio,clinical_variables)
        spath <- stabpath(y=y.tr, x=train_dataset, size = 0.632,
                          steps = 100,mc.cores=1, family="binomial",
                          weakness=.8, alpha=alph)
        stabproteinsidx <- c060::stabsel(spath, error=1, type="pfer", pi_thr=0.6)$stable
        names(stabproteinsidx) <- c(prot$AptName,"BMI","Age","Nulliparity")[stabproteinsidx]
        log.model <- glm(y.tr ~ .,  data=train_dataset[,stabproteinsidx], family = binomial)
        ratio <- x.te[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.te[,2505]
        clinical_variables <- x.te[,which(colnames(x.te)%in%c("BMI","parity","age"))]
        proteins <- x.te[,which(!colnames(x.te)%in%c("BMI","parity","age"))]
        test_dataset <- cbind(proteins,ratio,clinical_variables)
        ENSS.list.ratiomand[[nseed]][[k]] <- predict(log.model, newdata = test_dataset[,stabproteinsidx], type="response")
        ENSS.list.ratiomand[[nseed+length(nseeds)]][[k]] <- stabproteinsidx
        
      }
    }
    
    
    if(alph == 1){
      save(rf.list, file=paste0("~/data/rf",tp,"_15_151_34.RData"))
      save(logcl.list.ratio, file=paste0("~/data/logcl",tp,"ratio_15_151_34.RData"))
      save(logcl.list.bms, file=paste0("~/data/logcl",tp,"bms_15_151_34.RData"))
      save(logcl.list.ratio.bms, file=paste0("~/data/logcl",tp,"ratiobms_15_151_34.RData"))
      save(K.y,file=paste0("~/data/Ky",tp,"_15_151_34.RData"))
    }
    
    save(EN.list, file=paste0("~/data/EN",tp,"_alpha",alph,"_15_151_34.RData"))
    save(ENSS.list, file=paste0("~/data/ENSS",tp,"_alpha",alph,"_15_151_34.RData"))
    save(ENSS.list.ratio, file=paste0("~/data/ENSS",tp,"_alpha",alph,"ratio_15_151_34.RData"))
    save(ENSS.list.ratiomand, file=paste0("~/data/ENSS",tp,"_alpha",alph,"ratio_15_151_34.RData"))
    
    
  }
}














################## Find alphas for EN and ENSS ################
seedchosen <- 1
tps <- 1:3

alpha_EN <- c()
for ( tp in tps ){
  alphavec_EN <- c()
  for ( seed_i in seedchosen){
    alphaval_EN <- c()
    aucdf <- c()
    for ( alph in c(1,0.8,0.6,0.4,0.2) ){
      load(paste0("~/data/Ky",tp,"_15_151_34.RData"))
      load(paste0("~/data/EN",tp,"_alpha",alph,"_15_151_34.RData"))
      prediction_rocr <- prediction(EN.list[[seed_i]],K.y[[seed_i]])
      perf <- performance(prediction_rocr,"tpr","fpr")
      auc <- performance(prediction_rocr,"auc")
      aucdf <- cbind(aucdf,auc@y.values %>% unlist())
    }     
    
    colnames(aucdf) <- c(1,0.8,0.6,0.4,0.2)
    means <- apply(aucdf,2,mean)
    alphaval_EN <- c(alphaval_EN,which(means==max(means))[1] %>% names() %>% as.numeric())
    alphavec_EN <- cbind(alphavec_EN,alphaval_EN)
  }
  alpha_EN <- rbind(alpha_EN,alphavec_EN)
}


alpha_ENSS <- c()
for ( tp in 1:3){
  alphavec_ENSS <- c()
  for ( seed_i in seedchosen){ 
    alphaval_ENSS <- c()
    aucdf <- c()
    for ( alph in c(1,0.8,0.6,0.4,0.2) ){
      load(paste0("~/data/Ky",tp,"_15_151_34.RData"))
      load(paste0("~/data/ENSS",tp,"_alpha",alph,"_15_151_34.RData"))
      prediction_rocr <- prediction(ENSS.list[[seed_i]],K.y[[seed_i]])
      perf <- performance(prediction_rocr,"tpr","fpr")
      auc <- performance(prediction_rocr,"auc")
      aucdf <- cbind(aucdf,auc@y.values %>% unlist())
    }     
    
    colnames(aucdf) <- c(1,0.8,0.6,0.4,0.2)
    means <- apply(aucdf,2,mean)
    alphaval_ENSS <- c(alphaval_ENSS,which(means==max(means))[1] %>% names() %>% as.numeric())
    alphavec_ENSS <- cbind(alphavec_ENSS,alphaval_ENSS)
  }
  alpha_ENSS <- rbind(alpha_ENSS,alphavec_ENSS)
}






######### Plot prediction performances of EN/RF/ENSS models in Paper II ################
par(mfrow=c(2,3))

sub <- c("A. Visit 1","B. Visit 2", "C. Visit 3")
library(wesanderson)
colorsss <- wes_palette("Darjeeling1", 5, type = c("continuous"))
# par(mfrow=c(2,2))
for ( tp in 1:3 ){
  for ( seed_i in seedchosen){ 
    
    plot(c(),xlim=c(0,1),ylim=c(0,1),xlab="False positive rate",ylab="True positive rate")
    abline(a=0,b=1)
    mtext(sub[tp],side=3,line=1, at=-0.05)
    
    load(paste0("~/data/Ky",tp,"_15_151_34.RData"))
    
    load(paste0("~/data/rf",tp,"_15_151_34.RData"))
    plot(performance(prediction(lapply(rf.list[[seed_i]], function(x){x[,2]}),K.y[[seed_i]]),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col="royalblue3",spread.estimate="stderror")
    aucrf <- performance(prediction(lapply(rf.list[[seed_i]], function(x){x[,2]}),K.y[[seed_i]]),"auc")
    text(0.5,0.45,paste("mean AUC RF =",aucrf@y.values %>% unlist() %>% mean() %>% round(.,2)),col="royalblue3",adj=0)
    
    load(paste0("~/data/EN",tp,"_alpha",alpha_EN[tp,seed_i],"_15_151_34.RData"))
    plot(performance(prediction(EN.list[[seed_i]],K.y[[seed_i]]),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col=colorsss[3],spread.estimate="stderror")
    aucen <- performance(prediction(EN.list[[seed_i]],K.y[[seed_i]]),"auc")
    text(0.5,0.4,paste("mean AUC EN =",aucen@y.values %>% unlist() %>% mean() %>% round(.,2)),col=colorsss[3],adj=0)
    
    load(paste0("~/data/ENSS",tp,"_alpha",alpha_ENSS[tp,seed_i],"_15_151_34.RData"))
    plot(performance(prediction(ENSS.list[[seed_i]],K.y[[seed_i]]),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col=colorsss[1],spread.estimate="stderror")
    aucen <- performance(prediction(ENSS.list[[seed_i]],K.y[[seed_i]]),"auc")
    text(0.5,0.35,paste("mean AUC ENSS =",aucen@y.values %>% unlist() %>% mean() %>% round(.,2)),col=colorsss[1],adj=0)
    
    

    
  }
}








######### Comparisons between ratio and bms and ratio + bms in Paper II ##################

sub <- c("D. Visit 1","E. Visit 2", "F. Visit 3")
library(wesanderson)
colorsss <- wes_palette("Darjeeling1", 10, type = c("continuous"))
# par(mfrow=c(2,2))
for ( tp in 1:3 ){
  for ( seed_i in seedchosen){ 
    
    plot(c(),xlim=c(0,1),ylim=c(0,1),xlab="False positive rate",ylab="True positive rate")
    abline(a=0,b=1)
    mtext(sub[tp],side=3,line=1, at=-0.05)
    
    load(paste0("~/data/Ky",tp,"_15_151_34.RData"))
    
    load(paste0("~/data/logcl",tp,"ratio_15_151_34.RData"))
    plot(performance(prediction(logcl.list.ratio[[seed_i]],K.y[[seed_i]]),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col=colorsss[10],spread.estimate="stderror")
    auclog <- performance(prediction(logcl.list.ratio[[seed_i]],K.y[[seed_i]]),"auc")
    text(0.5,0.45,paste("mean AUC ratio + cl =",auclog@y.values %>% unlist() %>% mean() %>% round(.,2)),adj=0,col=colorsss[10])
    
    load(paste0("~/data/logcl",tp,"bms_15_151_34.RData"))
    plot(performance(prediction(logcl.list.bms[[seed_i]],K.y[[seed_i]]),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col=colorsss[1],spread.estimate="stderror")
    auclog <- performance(prediction(logcl.list.bms[[seed_i]],K.y[[seed_i]]),"auc")
    text(0.5,0.4,paste("mean AUC bms =",auclog@y.values %>% unlist() %>% mean() %>% round(.,2)),adj=0,col=colorsss[1])
    
    load(paste0("~/data/logcl",tp,"ratiobms_15_151_34.RData"))
    plot(performance(prediction(logcl.list.ratio.bms[[seed_i]],K.y[[seed_i]]),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col=6,spread.estimate="stderror")
    auclog <- performance(prediction(logcl.list.ratio.bms[[seed_i]],K.y[[seed_i]]),"auc")
    text(0.5,0.35,paste("mean AUC ratio + bms =",auclog@y.values %>% unlist() %>% mean() %>% round(.,2)),adj=0,col=6)
    
  }
}







############### Performance of ENSS only for Paper II #####################

sub <- c("Visit 1","Visit 2", "Visit 3")
library(wesanderson)
colorsss <- wes_palette("Darjeeling1", 5, type = c("continuous"))
par(mfrow=c(2,2))
for ( tp in 1:3 ){
  for ( seed_i in seedchosen){ 
    
    plot(c(),xlim=c(0,1),ylim=c(0,1),xlab="False positive rate",ylab="True positive rate")
    abline(a=0,b=1)
    mtext(sub[tp],side=3,line=1, at=-0.05)
    
    load(paste0("~/data/Ky",tp,"_15_151_34.RData"))
    
    load(paste0("~/data/ENSS",tp,"_alpha",alpha_ENSS[tp,seed_i],"_15_151_34.RData"))
    plot(performance(prediction(ENSS.list[[seed_i]],K.y[[seed_i]]),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col="orange",spread.estimate="boxplot")
    aucen <- performance(prediction(ENSS.list[[seed_i]],K.y[[seed_i]]),"auc")
    text(0.6,0.35,paste("AUC ENSS =",aucen@y.values %>% unlist() %>% mean() %>% round(.,2)),col="orange",adj=0)
    
    
  }
}






















