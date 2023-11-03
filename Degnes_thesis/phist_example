###########
#
# Author: Maren-Helene Langeland Degnes 
# Date: 27. Nov 2022
# Purpose: Plot p-value histograms examples for my thesis
#
###########


par(mfrow=c(1,3))
set.seed(1)
hist(sample(seq(from=0,to=1,length.out=1000),5000,replace=T),
     breaks=seq(0,1,0.05),ylim=c(0,800),
     main="A. Keep the global null hypothesis",
     xlab="P-values")
text(0.5,850,"No true positive results",xpd=T)
text(0.5,800,"Assumptions are met",xpd=T)
abline(h=qbinom(0.975,5000,0.05),col="red")            # Code from Breheny et al. (2018)
abline(h=qbinom(1-0.05*0.05,5000,0.05),col="blue")     # Code from Breheny et al. (2018)

hist(sample(c(seq(from=0,to=1,length.out=1000),
              seq(from=0.6,to=0.7,length.out=70),
              seq(from=0.55,to=0.67,length.out=50),
              seq(from=0.7,to=0.75,length.out=40)),5000,replace=T),
     breaks=seq(0,1,0.05),ylim=c(0,800),
     main="B. Keep the global null hypothesis",
     xlab="P-values")
text(0.5,850,"No true positive results",xpd=T)
text(0.5,800,"Assumptions are not met",xpd=T)
abline(h=qbinom(0.975,5000,0.05),col="red")            # Code from Breheny et al. (2018)
abline(h=qbinom(1-0.05*0.05,5000,0.05),col="blue")     # Code from Breheny et al. (2018)

hist(sample(c(seq(from=0,to=1,length.out=1000),
              seq(from=0,to=0.05,length.out=100)),5000,replace=T),
     breaks=seq(0,1,0.05),ylim=c(0,800),
     main="C. Reject the global null hypothesis",
     xlab="P-values")
text(0.5,850,"True positive results are likely",xpd=T)
text(0.5,800,"Assumptions are met",xpd=T)
abline(h=qbinom(0.975,5000,0.05),col="red")            # Code from Breheny et al. (2018)
abline(h=qbinom(1-0.05*0.05,5000,0.05),col="blue")     # Code from Breheny et al. (2018)

# References: Breheny, P., A. Stromberg, and J. Lambert. 2018. 'p-Value Histograms: Inference and Diagnostics', High Throughput, 7.
