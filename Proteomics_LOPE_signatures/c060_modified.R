### 
#
# code modified by Maren-Helene Langeland Degnes February 2023
# Purpose of modification is to let clinical variables always be included in the elastic net with stability selection models
#
###


stabpath <- function(y,x,size=0.632,steps=100,weakness=1,mc.cores=getOption("mc.cores", 2L),...){
  
  require(c060)
  require(glmnet)
  require(parallel)
  
  fit <- glmnet(x,y,...)
  if(class(fit)[1]=="multnet"|class(fit)[1]=="lognet") y <- as.factor(y)
  #if(class(fit)[1]=="lognet") y <- as.logical(y) 
  p <- ncol(x)
  #draw subsets
  subsets <- sapply(1:steps,function(v){sample(1:nrow(x),nrow(x)*size)})
  
  # parallel computing depending on OS
  # UNIX/Mac
  if (mc.cores > 1){
    if (.Platform$OS.type != "windows") {
      res <- mclapply(1:steps, mc.cores = mc.cores, glmnet.subset, 
                      subsets, x, y, lambda = fit$lambda, weakness, p, 
                      ...)
    }
    else {
      # Windows  
      cl  <- makePSOCKcluster(mc.cores)
      clusterExport(cl,c("glmnet","drop0"))
      res <- parLapply(cl, 1:steps,glmnet.subset,subsets,
                       x,y,lambda=fit$lambda,weakness,p,...)
      stopCluster(cl)
    }
  }else{
    print("test: mc.cores=1")
    res <- lapply(1:steps, glmnet.subset, 
                  subsets, x, y, lambda = fit$lambda, weakness, p, ...)
  }
  
  #merging
  res <- res[unlist(lapply(lapply(res,dim),function(x) x[2]==dim(res[[1]])[2]))]
  x <- as.matrix(res[[1]])
  qmat <- matrix(ncol=ncol(res[[1]]),nrow=length(res))
  qmat[1,] <- colSums(as.matrix(res[[1]]))
  for(i in 2:length(res)){
    qmat[i,] <- colSums(as.matrix(res[[i]]))
    x <- x + as.matrix(res[[i]])
  }
  x <- x/length(res)
  qs <- colMeans(qmat)
  out <- list(fit=fit,x=x,qs=qs)	
  class(out) <- "stabpath" 
  return(out)
}


#internal function used by lapply 
glmnet.subset <- function(index,subsets,x,y,lambda,weakness,p,...){
  if(length(dim(y))==2|class(y)=="Surv"){
    glmnet(x[subsets[,index],],y[subsets[,index],],lambda=lambda
           ,penalty.factor= c(1/runif(p-4,weakness,1),rep(0,4)),...)$beta!=0              ##### I changed the penalty.factor argument from
  }else{                                                                                  ##### 1/runif(p,weaknes,1) to this inside c()
    if(is.factor(y)&length(levels(y))>2){
      temp <- glmnet(x[subsets[,index],],y[subsets[,index]],lambda=lambda
                     ,penalty.factor= c(1/runif(p-4,weakness,1),rep(0,4)),...)[[2]]       ##### Same change of penalty factor as above
      temp <- lapply(temp,as.matrix)
      Reduce("+",lapply(temp,function(x) x!=0))>0
      
    }	
    else{
      glmnet(x[subsets[,index],],y[subsets[,index]],lambda=lambda
             ,penalty.factor= c(1/runif(p-4,weakness,1),rep(0,4)),...)$beta!=0            ##### Same change of penalty factor as above
    }
  }	
}

