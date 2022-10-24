#########
#
# Author: Maren-Helene Langeland Degnes
# Date: 24.Oct.22
# 
#
########




library(tidyverse)
library(pheatmap)
library(readxl)
library(gridExtra)
library(ggrepel)



par(mfrow=c(2,2),mar=c(4,4,4,4),cex.lab=1.5,cex.main=1.5)
p <- rep(NA, length(MatNormTest$`p value`))   # vector with same length as p-values/proteins
b <- 0.05   # bin width



load("placentareltak_testres.RData")
hist(MatNormTest$`p value`, 
     breaks=seq(0,1,b), 
     main=paste0("A. 4V Placental uptake and release"),
     xlab="P-value")
legend("topright",legend=paste("q<0.05:",length(which(MatNormTest$`adjusted p`<0.05)),"proteins"),cex=1.2)

# Higher criticism and quality control thresholds - code from Breheny et al. 2018 (Reference in bottom of script)
abline(h=qbinom(0.975,length(MatNormTest$`p value`),b),col="red")
abline(h=qbinom(1-b*0.05,length(MatNormTest$`p value`),b),col="blue")




T1 <- read_excel("T1testres.xlsx")
T2 <- read_excel("T2testres.xlsx")
T3 <- read_excel("T3testres.xlsx")

hist(T1$P.Value,breaks = seq(0,1,by=b),main="B. Visit 1",xlab="P-value")
abline(h=qbinom(0.975,length(p),b),col="red")
abline(h=qbinom(1-b*0.05,length(p),b),col="blue")
legend("topright",legend=paste("q<0.05:",length(which(T1$adj.P.Val<0.05)),"proteins"),cex=1.2)

hist(T2$P.Value,breaks = seq(0,1,by=b),main="C. Visit 2",xlab="P-value",ylim=c(0,320))
abline(h=qbinom(0.975,length(p),b),col="red")
abline(h=qbinom(1-b*0.05,length(p),b),col="blue")
legend("topright",legend=paste("q<0.05:",length(which(T2$adj.P.Val<0.05)),"proteins"),cex=1.2)

hist(T3$P.Value,breaks = seq(0,1,by=b),main="D. Visit 3",xlab="P-value")
abline(h=qbinom(0.975,length(p),b),col="red")
abline(h=qbinom(1-b*0.05,length(p),b),col="blue")
legend("topright",legend=paste("q<0.05:",length(which(T3$adj.P.Val<0.05)),"proteins"),cex=1.2)



Breheny, P., A. Stromberg, and J. Lambert. 2018. 'p-Value Histograms: Inference and Diagnostics', High Throughput, 7.
