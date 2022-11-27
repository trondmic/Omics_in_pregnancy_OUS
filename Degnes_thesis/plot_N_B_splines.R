#######
#
# Author: Maren-Helene Langeland Degnes
# Date: 24. Oct. 2022
# Purpose: Plot natural and B splines functions
#
#########
library(tidyverse)

GA <- c(rnorm(70,mean=15.5,sd=1),
        rnorm(70,mean=24,sd=1),
        rnorm(70,mean=33,sd=0.7)) %>% sort()


par(mfrow=c(1,2))
## Plotting GA against natural splines functions of GA
plot(c(),xlim=range(GA),ylim=c(-0.5,1.5),                            # empty plot
     main="A. Natural splines functions - linear at boundary knots",
     ylab="f(Gestational weeks)",xlab="Gestational weeks")
naturalsplinesmx <- ns(GA,                                           # create natural splines functions
                       knots=15.5,
                       Boundary.knots=c(15,34))
matlines(GA, naturalsplinesmx)			     # fill plot with GA vs natural splines functions
lines(GA, rowSums(naturalsplinesmx), col = "purple",lwd=2)           # calculate sum of splines and plot
abline(v=c(15,15.5,34),col="cornflowerblue",lty=2)                           # lines for knots

## Plotting GA against B splines functions of GA
plot(c(),xlim=range(GA),ylim=c(-0.5,1.5),                           # empty plot
     main="B. B splines functions - sum to 1 inside inner knots",
     ylab="f(Gestational weeks)",xlab="Gestational weeks")
bsplinesmx <- bs(GA,                                                # create B splines functions
                 degree=2,
                 knots=15.5,
                 Boundary.knots=c(15,34))
matlines(GA, bsplinesmx) 		         	    # fill plot with GA vs natural splines functions
lines(GA, rowSums(bsplinesmx), col = "purple",lwd=2)                 # calculate sum of splines and plot 
abline(v=c(15,15.5,34),col="cornflowerblue",lty=2)                           # lines for knots
