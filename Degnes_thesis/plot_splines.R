#######
#
# Author: Maren-Helene Langeland Degnes
# Date: 24.Oct.22
#
#########

#GA <- read.xlsx("GAs.xlsx")

par(mfrow=c(1,2))
## Plotting GA against natural splines functions of GA
plot(c(),xlim=range(GA),ylim=c(-1.5,1.5),                            # empty plot
     main="Natural splines functions - linear at boundary knots",
     ylab="f(GA)",xlab="GA")
naturalsplinesmx <- ns(GA,                                           # create natural splines functions
                       knots=15.5,
                       Boundary.knots=c(15,34))
matlines(GA, naturalsplinesmx)			     # fill plot with GA vs natural splines functions
lines(GA, rowSums(naturalsplinesmx), col = "purple",lwd=2)           # calculate sum of splines and plot
abline(v=c(15,15.5,34),col="green",lty=2)                           # lines for knots

## Plotting GA against B splines functions of GA
plot(c(),xlim=range(GA),ylim=c(-1.5,1.5),                           # empty plot
     main="B splines functions - sum to 1 inside inner knots",
     ylab="f(GA)",xlab="GA")
bsplinesmx <- bs(GA,                                                # create B splines functions
                 degree=2,
                 knots=15.5,
                 Boundary.knots=c(15,34))
matlines(GA, bsplinesmx) 		         	    # fill plot with GA vs natural splines functions
lines(GA, rowSums(bsplinesmx), col = "purple",lwd=2)                 # calculate sum of splines and plot abline(v=c(15,15.5,34),col="green",lty=2)                           # lines for knots
