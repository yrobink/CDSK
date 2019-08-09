
#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

###############
## Libraries ##
###############

library(devtools)
load_all( "../R/CDSK" )
library(rgl)

##########
## main ##
##########


l63 = CDSK::Lorenz63$new( size = 10 )
X0  = l63$orbit( base::seq( 0 , 10 , length = 100 ) )
X   = X0[,,10]

plot3d( X[,1] , X[,2] , X[,3] , col = "red" )

base::cat("Done\n")
