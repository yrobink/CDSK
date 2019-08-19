
#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

base::rm( list = base::ls() )

###############
## Libraries ##
###############

library(devtools)
load_all( "../R/CDSK" )
library(rgl)

##########
## main ##
##########


l63 = CDSK::Lorenz63$new( size = 1000 )
X0  = l63$orbit( base::seq( 0 , 400 , length = 1000 ) )
X   = X0[1000,,]

plot3d( X[,1] , X[,2] , X[,3] , col = "red" )

base::cat("Done\n")
