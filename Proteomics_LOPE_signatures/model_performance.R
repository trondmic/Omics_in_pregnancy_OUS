#####
# 
# Author: Maren-Helene Langeland Degnes
# Date: feb-2023
#
# Content: Predict LOPE using random forest, elastic net and logistic 
#          regression algorithms results can be found in Figure 2, 
#          elastic net with stability selection is also included in the for-loop
#
#     A. training and testing in 5 fold CV
#     B. Find alpha that gave the highest performing elastic net
#     C. Plot performances of all models into Figure 3
#
#     I moderated the functions from the c060 package so the clinical variables
#     could be included in every elastic net model as mandatory variables. 
#     The moderated code can be found here: 
#     https://github.com/trondmic/Omics_in_pregnancy_OUS/blob/main/Proteomics_LOPE_signatures/c060_modified.R
#
#
#####




# Load packages
library(tidyverse)
library(randomForest)
library(ROCR)
library(glmnet)
library(Biobase)
library(parallel)
library(c060)
source("scripts/04_c060_functions.R")  # Load c060 functions that I moderated





############### A. Train and test a diversity of machine learning models in a 5-fold cross validation ##################

## Load data and sample train and test ####
load("data/MoM_STORK_log2win_exprset.RData")
MoMmx <- exprs(MOM_STORK_log2win_exprset) %>% t()
ano <- pData(MOM_STORK_log2win_exprset)
prot <- fData(MOM_STORK_log2win_exprset)
colnames(MoMmx) <- prot$AptName


MoMmx <- MoMmx[which(ano$PEgroups %in% c(0,2)),]  # Select participants
ano <- ano[which(ano$PEgroups %in% c(0,2)),]
ranges <- tapply(ano$GAWeeks , ano$TimePoint, range) # save ranges
table(rownames(MoMmx) == rownames(ano))   # Check that matrices dimensions and names match
ano$FOLDID <- NA  # Create variable FOLDID in which I separate participants into 5 folds for cross validation
all_ids <- ano$ID %>% unique()   # Create vector with all participant ids
set.seed(1,kind="default")                            # set seed so the sampling is reproducible 
for ( id in all_ids ){  ano$FOLDID[which(ano$ID == id)] <- sample(1:5,1)} # separate into folds
table(ano$FOLDID,ano$PEgroups)                        # print the distribution of participant groups in each fold


tps <- 1:3
FOLDS <- 1 # We ran the whole analysis multiple times 
for ( tp in 1:3){  # For each time point, fit all the different models 
  
  tp.latePE.idx <- which(ano$TimePoint==tp)  # Choose data on time point tp
  y_mx <- ano[tp.latePE.idx,]
  x_mx <- MoMmx[tp.latePE.idx,]
  
  table(rownames(x_mx) == rownames(y_mx))  # Check that dimensions and names match between annotation and features matrices
  
  x <- data.frame(x_mx,BMI=y_mx$BMI,age=as.numeric(y_mx$Age),parity=y_mx$Nulliparity%>%as.factor())   # Set clinical data and number of folds for CV
  x_nocl <- data.frame(x_mx)
  y <- y_mx$PE %>% as.factor()
  n <- length(y)
  K <- 5
  
  for ( alph in c(1,0.8,0.6,0.4,0.2) ){ # For each value of alpha, fit an elastic net model 
    
    rf.list <- list(); for (i in FOLDS){rf.list[[i]] <- list()}
    EN.list <- list(); for (i in c(FOLDS,FOLDS+length(FOLDS))){EN.list[[i]] <- list()}
    ENSS.list <- list(); for (i in c(FOLDS,FOLDS+length(FOLDS))){ENSS.list[[i]] <- list()}
    ENSS.list.ratio <- list(); for (i in c(FOLDS,FOLDS+length(FOLDS))){ENSS.list.ratio[[i]] <- list()}
    ENSS.list.ratiomand <- list(); for (i in c(FOLDS,FOLDS+length(FOLDS))){ENSS.list.ratiomand[[i]] <- list()}
    logcl.list.ratio <- list(); for (i in FOLDS){logcl.list.ratio[[i]] <- list()}
    logcl.list.bms <- list(); for (i in FOLDS){logcl.list.bms[[i]] <- list()}
    logcl.list.ratio.bms <- list(); for (i in FOLDS){logcl.list.ratio.bms[[i]] <- list()}
    K.y <- list(); for (i in FOLDS){K.y[[i]] <- list()}
    
    for (FOLDID in FOLDS){
      
      # Run training and testing
      for (k in 1:K){
        
        x.tr <- x[-which(y_mx[,47+FOLDID] == k),] %>% as.data.frame()
        x.te <- x[which(y_mx[,47+FOLDID] == k),] %>% as.data.frame()
        
        x.tr.nocl <- x_nocl[-which(y_mx[,47+FOLDID] == k),] %>% as.data.frame()
        x.te.nocl <- x_nocl[which(y_mx[,47+FOLDID] == k),] %>% as.data.frame()
        
        y.tr <- y[-which(y_mx[,47+FOLDID] == k)] %>% as.factor()
        y.te <- y[which(y_mx[,47+FOLDID] == k)] %>% as.factor()
        K.y[[FOLDID]][[k]] <- y.te
        
        
        
        if(alph == 1){          #Run random forests only once, not for every alpha
          
          # i Random forests!
          set.seed(500)
          rf.mod <- randomForest(x= x.tr,  y= y.tr,  ntree=1000, keep.forest=TRUE, importance=TRUE)
          rf.list[[FOLDID]][[k]] <- predict(rf.mod,  newdata= x.te,  type="prob")
          
          
          # ii Logistic regression including FLT1/PGF ratio + clinical variables
          # Results are found in Figure 3D-E, blue line
          ratio <- x.tr[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.tr[,2505]
          clinical_variables <- x.tr[,which(colnames(x.tr)%in%c("BMI","parity","age"))]
          train_dataset <- cbind(ratio,clinical_variables)
          log.model <- glm(y.tr ~ .,  data=train_dataset, family = binomial)
          ratio <- x.te[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.te[,2505]
          clinical_variables <- x.te[,which(colnames(x.te)%in%c("BMI","parity","age"))]
          test_dataset <- cbind(ratio, clinical_variables)
          logcl.list.ratio[[FOLDID]][[k]] <- predict(log.model, newdata = test_dataset, type="response")
          
          
          # iii Logistic regression including the biomarker signature
          # Results are found in Figure 3D-E, red line
          clinical_variables <- x.tr[,which(colnames(x.tr)%in%c("BMI","parity","age"))]
          biomarker_signature_proteins <- x.tr[,list(c(),c(1319,3159,3797),c(1480, 1544, 3159, 3928, 4008))[[tp]]]
          train_dataset <- cbind(biomarker_signature_proteins,clinical_variables)
          log.model <- glm(y.tr ~ .,  data=train_dataset, family = binomial)
          clinical_variables <- x.te[,which(colnames(x.te)%in%c("BMI","parity","age"))]
          biomarker_signature_proteins <- x.te[,list(c(),c(1319,3159,3797),c(1480, 1544, 3159, 3928, 4008))[[tp]]]
          test_dataset <- cbind(biomarker_signature_proteins, clinical_variables)
          logcl.list.bms[[FOLDID]][[k]] <- predict(log.model, newdata = test_dataset, type="response")
          
          
          # iv Logistic regression including the biomarker signature + FLT1/PGF ratio + clinical variables
          # Results are found in Figure 3D-E, purple line
          ratio <- x.tr[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.tr[,2505]
          clinical_variables <- x.tr[,which(colnames(x.tr)%in%c("BMI","parity","age"))]
          biomarker_signature_proteins <- x.tr[,list(c(),c(1319,3159,3797),c(1480, 1544, 3159, 3928, 4008))[[tp]]]
          train_dataset <- cbind(ratio,biomarker_signature_proteins,clinical_variables)
          log.model <- glm(y.tr ~ .,  data=train_dataset, family = binomial)
          ratio <- x.te[,which(prot$EntrezGeneSymbol=="FLT1")[2]] - x.te[,2505]
          clinical_variables <- x.te[,which(colnames(x.te)%in%c("BMI","parity","age"))]
          biomarker_signature_proteins <- x.te[,list(c(),c(1319,3159,3797),c(1480, 1544, 3159, 3928, 4008))[[tp]]]
          test_dataset <- cbind(ratio, biomarker_signature_proteins, clinical_variables)
          logcl.list.ratio.bms[[FOLDID]][[k]] <- predict(log.model, newdata = test_dataset, type="response")
          
        }
        
        set.seed(500)
        # v Run elastic net, including LASSO (alpha=1)
        EN.mod <- cv.glmnet(x=x.tr %>% data.matrix(),y=y.tr,alpha=alph,type.measure="class",
                            family="binomial",penalty.factor=c(rep(1,ncol(MoMmx)),rep(0,3)))
        EN.list[[FOLDID]][[k]] <- predict(EN.mod$glmnet.fit,s=EN.mod$lambda.min,
                                          newx=x.te%>%data.matrix(), type="response")
        EN.coef <- predict(EN.mod, s = EN.mod$lambda.min , type = "coefficients")
        EN.coef <- as.matrix(EN.coef)                       # Prepare the prediction for matrix-indexing
        EN.list[[FOLDID+length(FOLDS)]][[k]] <- EN.coef[EN.coef != 0,]
        
        
        # vi Run elastic net with stability selection and then logistic regression
        set.seed(500)
        source("scripts/04_c060_functions.R")
        spath <- stabpath(y=y.tr, x=x.tr, size = 0.632,
                          steps = 100,mc.cores=1, family="binomial",
                          weakness=.8, alpha=alph)
        stabproteinsidx <- c060::stabsel(spath, error=1, type="pfer", pi_thr=0.6)$stable
        names(stabproteinsidx) <- c(prot$AptName,"BMI","Age","Nulliparity")[stabproteinsidx]
        log.model <- glm(y.tr ~ .,  data=x.tr[,stabproteinsidx], family = binomial)
        ENSS.list[[FOLDID]][[k]] <- predict(log.model, newdata = x.te[,stabproteinsidx], type="response")
        ENSS.list[[FOLDID+length(FOLDS)]][[k]] <- stabproteinsidx
        
      }
    }
    
    # Save all results objects!
    if(alph == 1){
      save(rf.list, file=paste0("data/rf",tp,".RData"))
      save(logcl.list.ratio, file=paste0("data/logcl",tp,"ratio.RData"))
      save(logcl.list.bms, file=paste0("data/logcl",tp,"bms.RData"))
      save(logcl.list.ratio.bms, file=paste0("data/logcl",tp,"ratiobms.RData"))
      save(K.y,file=paste0("data/Ky",tp,".RData")) # true groups vector
    }
    
    save(EN.list, file=paste0("data/EN",tp,"_alpha",alph,".RData"))
    save(ENSS.list, file=paste0("data/ENSS",tp,"_alpha",alph,".RData"))
    
    
  }
}












############# B. Find alpha that gave the highest performing elastic net per time point ##

FOLDchosen <- 1
tps <- 1:3

# Elastic net
alpha_EN <- c()
for ( tp in tps ){
  alphavec_EN <- c()
  for ( FOLD in FOLDchosen){
    alphaval_EN <- c()
    aucdf <- c()
    for ( alph in c(1,0.8,0.6,0.4,0.2) ){
      load(paste0("data/Ky",tp,".RData"))
      load(paste0("data/EN",tp,"_alpha",alph,".RData"))
      prediction_rocr <- prediction(EN.list[[FOLD]],K.y[[FOLD]])
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

# Elastic net with stability selection
alpha_ENSS <- c()
for ( tp in 1:3){
  alphavec_ENSS <- c()
  for ( FOLD in FOLDchosen){ 
    alphaval_ENSS <- c()
    aucdf <- c()
    for ( alph in c(1,0.8,0.6,0.4,0.2) ){
      load(paste0("data/Ky",tp,".RData"))
      load(paste0("data/ENSS",tp,"_alpha",alph,".RData"))
      prediction_rocr <- prediction(ENSS.list[[FOLD]],K.y[[FOLD]])
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










######### Figure 3 ################

jpeg("Figure3.jpg",width=3750,height=2500,res=300)
par(mfrow=c(2,3))
minmax_tab <- matrix(nrow=3,ncol=3)
CI_tab <- matrix(nrow=6,ncol=3)
colnames(minmax_tab) <- c("Tp1","Tp2","Tp3")

sub <- c("a. Visit 1","b. Visit 2", "c. Visit 3")
library(wesanderson)
colorsss <- wes_palette("Darjeeling1", 5, type = c("continuous"))
# par(mfrow=c(2,2))
for ( tp in 1:3 ){
  for ( FOLD in FOLDchosen){ 
    
    plot(c(),xlim=c(0,1),ylim=c(0,1),xlab="False positive rate",ylab="True positive rate")
    abline(a=0,b=1)
    mtext(sub[tp],side=3,line=1, at=-0.05)
    
    load(paste0("data/Ky",tp,".RData"))
    y <- K.y[[FOLD]]
    
    load(paste0("data/rf",tp,".RData"))
    class_prob <- lapply(rf.list[[FOLD]], function(x){x[,2]})
    plot(performance(prediction(class_prob,y),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col="royalblue3",spread.estimate="stderror")
    aucrf <- performance(prediction(class_prob,y),"auc")
    minmax_tab[1,tp] <- paste(aucrf@y.values %>% unlist() %>% mean() %>% round(.,2),aucrf@y.values %>% unlist() %>% min() %>% round(.,2), aucrf@y.values %>% unlist() %>% max() %>% round(.,2))
    CI_vec <- ci.auc(roc(unlist(y),unlist(class_prob)))
    CI_tab[1,tp] <- paste0(aucrf@y.values %>% unlist() %>% mean() %>% round(.,2), " (", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")
    text(0.25,0.2,paste("AUC RF =",aucrf@y.values %>% unlist() %>% mean() %>% round(.,2),paste0("(", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")),col="royalblue3",adj=0)
    
    load(paste0("data/EN",tp,"_alpha",alpha_EN[tp,FOLD],".RData"))
    class_prob <- EN.list[[FOLD]]
    plot(performance(prediction(class_prob,y),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,
         col="darkorange1",spread.estimate="stderror")
    aucen <- performance(prediction(class_prob,y),"auc")
    minmax_tab[2,tp] <- paste(aucen@y.values %>% unlist() %>% mean() %>% round(.,2), aucen@y.values %>% unlist() %>% min() %>% round(.,2), aucen@y.values %>% unlist() %>% max() %>% round(.,2))
    CI_vec <- ci.auc(roc(unlist(y),unlist(class_prob)))
    CI_tab[2,tp] <- paste0(aucen@y.values %>% unlist() %>% mean() %>% round(.,2), " (", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")
    text(0.25,0.15,paste("AUC EN =",aucen@y.values %>% unlist() %>% mean() %>% round(.,2),paste0("(", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")),
         col="darkorange1",adj=0)
    
    load(paste0("data/ENSS",tp,"_alpha",alpha_ENSS[tp,FOLD],".RData"))
    class_prob <- ENSS.list[[FOLD]]
    plot(performance(prediction(class_prob,y),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col=colorsss[1],spread.estimate="stderror")
    aucen <- performance(prediction(class_prob,y),"auc")
    minmax_tab[3,tp] <-  paste(aucen@y.values %>% unlist() %>% mean() %>% round(.,2), aucen@y.values %>% unlist() %>% min() %>% round(.,2), aucen@y.values %>% unlist() %>% max() %>% round(.,2))
    CI_vec <- ci.auc(roc(unlist(y),unlist(class_prob)))
    CI_tab[3,tp] <- paste0(aucen@y.values %>% unlist() %>% mean() %>% round(.,2), " (", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")
    text(0.25,0.1,paste("AUC ENSS =",aucen@y.values %>% unlist() %>% mean() %>% round(.,2),paste0("(", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")),col=colorsss[1],adj=0)
    
    
  }
}





######### Comparisons between ratio and bms and ratio + bms ##################

sub <- c("d. Visit 1","e. Visit 2", "f. Visit 3")
minmax_tab <- matrix(nrow=3,ncol=3)
colnames(minmax_tab) <- c("Tp1","Tp2","Tp3")

library(wesanderson)
colorsss <- wes_palette("Darjeeling1", 10, type = c("continuous"))
# par(mfrow=c(2,2))
for ( tp in 1:3 ){
  for ( FOLD in FOLDchosen){ 
    
    plot(c(),xlim=c(0,1),ylim=c(0,1),xlab="False positive rate",ylab="True positive rate")
    abline(a=0,b=1)
    mtext(sub[tp],side=3,line=1, at=-0.05)
    
    load(paste0("data/Ky",tp,".RData"))
    y <- K.y[[FOLD]]
    
    load(paste0("data/logcl",tp,"ratio.RData"))
    class_prob <- logcl.list.ratio[[FOLD]]
    plot(performance(prediction(class_prob,y),"tpr","fpr"),lwd=3,
         avg="vertical",add=TRUE,col="dodgerblue",spread.estimate="stderror")
    auclog <- performance(prediction(class_prob,y),"auc")
    minmax_tab[1,tp] <-  paste(auclog@y.values %>% unlist() %>% mean() %>% round(.,2), auclog@y.values %>% unlist() %>% min() %>% round(.,2), auclog@y.values %>% unlist() %>% max() %>% round(.,2))
    CI_vec <- ci.auc(roc(unlist(y),unlist(class_prob)))
    CI_tab[4,tp] <- paste0(auclog@y.values %>% unlist() %>% mean() %>% round(.,2), " (", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")
    text(0.25,0.2,paste("AUC ratio + cl =",auclog@y.values %>% unlist() %>% mean() %>% round(.,2),paste0("(", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")),adj=0,
         col="dodgerblue")
    
    load(paste0("data/logcl",tp,"bms.RData"))
    class_prob <- logcl.list.bms[[FOLD]]
    plot(performance(prediction(class_prob,y),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col=colorsss[1],spread.estimate="stderror")
    auclog <- performance(prediction(class_prob,y),"auc")
    minmax_tab[2,tp] <-  paste(auclog@y.values %>% unlist() %>% mean() %>% round(.,2), auclog@y.values %>% unlist() %>% min() %>% round(.,2), auclog@y.values %>% unlist() %>% max() %>% round(.,2))
    CI_vec <- ci.auc(roc(unlist(y),unlist(class_prob)))
    CI_tab[5,tp] <- paste0(auclog@y.values %>% unlist() %>% mean() %>% round(.,2), " (", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")
    text(0.25,0.15,paste("AUC biomarker signature =",auclog@y.values %>% unlist() %>% mean() %>% round(.,2),paste0("(", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")),adj=0,col=colorsss[1])
    
    load(paste0("data/logcl",tp,"ratiobms.RData"))
    class_prob <- logcl.list.ratio.bms[[FOLD]]
    plot(performance(prediction(class_prob,y),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col=6,spread.estimate="stderror")
    auclog <- performance(prediction(class_prob,y),"auc")
    minmax_tab[3,tp] <-  paste(auclog@y.values %>% unlist() %>% mean() %>% round(.,2), auclog@y.values %>% unlist() %>% min() %>% round(.,2), auclog@y.values %>% unlist() %>% max() %>% round(.,2))
    CI_vec <- ci.auc(roc(unlist(y),unlist(class_prob)))
    CI_tab[6,tp] <- paste0(auclog@y.values %>% unlist() %>% mean() %>% round(.,2), " (", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")
    text(0.25,0.1,paste("AUC ratio + biomarker signature =",auclog@y.values %>% unlist() %>% mean() %>% round(.,2),paste0("(", CI_vec[1] %>% round(.,2), "-", CI_vec[3] %>% round(.,2),")")),adj=0,col=6)
    
  }
}
dev.off()


write.csv2(CI_tab,file="ScientRep/CI.csv")









############### Comparisons between RF, EN and ENSS #####################

sub <- c("Visit 1","Visit 2", "Visit 3")
library(wesanderson)
colorsss <- wes_palette("Darjeeling1", 5, type = c("continuous"))
par(mfrow=c(2,2))
for ( tp in 1:3 ){
  for ( FOLD in FOLDchosen){ 
    
    plot(c(),xlim=c(0,1),ylim=c(0,1),xlab="False positive rate",ylab="True positive rate")
    abline(a=0,b=1)
    mtext(sub[tp],side=3,line=1, at=-0.05)
    
    load(paste0("data/Ky",tp,".RData"))
    
    load(paste0("data/ENSS",tp,"_alpha",alpha_ENSS[tp,FOLD],".RData"))
    plot(performance(prediction(ENSS.list[[FOLD]],K.y[[FOLD]]),"tpr","fpr"),lwd=3,avg="vertical",add=TRUE,col="orange",spread.estimate="boxplot")
    aucen <- performance(prediction(ENSS.list[[FOLD]],K.y[[FOLD]]),"auc")
    text(0.6,0.35,paste("AUC ENSS =",aucen@y.values %>% unlist() %>% mean() %>% round(.,2)),col="orange",adj=0)
    
    
  }
}

