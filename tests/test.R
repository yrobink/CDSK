
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
library(plot3D)


###############
## Functions ##
###############

## How SDCK run {{{

lorenz63 = function( t , X , parms = list( s = 10 , r = 28 , b = 2.667 , size = 1 ) )
{
	s = parms$s
	r = parms$r
	b = parms$b
	size = parms$size
	dim  = 3
	i = list()
	for( k in 1:dim )
	{
		i[[k]] = base::seq( k , dim * size , dim )
	}
	
	dX = numeric( length(X) )
	dX[i[[1]]] = s * ( X[i[[2]]] - X[i[[1]]] )
	dX[i[[2]]] = r * X[i[[1]]] - X[i[[2]]] - X[i[[1]]] * X[i[[3]]]
	dX[i[[3]]] = X[i[[1]]] * X[i[[2]]] - b * X[i[[3]]]
	return( list(dX) )
}

random_IC = function( box_ic , n_ic = 1 )
{
	dim = base::dim(box_ic)[1]
	ic  = stats::runif( dim * n_ic , min = box_ic[,1] , box_ic[,2] )
	return(ic)
}

build_lorenz63 = function( t , X0 = NULL , size = 1 , s = 10 , r = 28 , b = 2.667 )
{
	box_ic = base::matrix( base::c( -20 , -20 , 0 , 20 , 20 , 40 ) , nrow = 3 , ncol = 2 )
	dim  = 3
	
	X0 = if( !is.null(X0) ) X0 else random_IC( box_ic , size )
	
	X  = deSolve::ode( X0 , t , lorenz63 , parms = list( s = s , r = r , b = b , size = size ) , method = "rk4" )[,-1]
	
	if( size > 1 )
	{
		X = base::aperm( base::array( X , dim = base::c( length(t) , dim , size ) ) , perm = base::c( 1 , 3 , 2 ) )
	}
	return(X)
}

##}}}

test_lorenz63 = function( plot = FALSE )##{{{
{
	t  = base:::seq( 0 , 50 , length = 1000 )
	size = 1000
	l63 = Lorenz63$new( size = size )
	X   = l63$orbit( t = t )
	
	
	if( plot )
	{
		grDevices::dev.new()
		graphics::par( mfrow = base::c( 2 , 2 ) )
		
		## Plot0
		snap = 1
		orb  = 1:size
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "red" , main = "Initial snapshot" , ticktype = "detailed" )
		
		## Plot1
		snap = length(t)
		orb  = 1:size
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "red" , main = "Final snapshot" , ticktype = "detailed" )
		
		## Plot3
		snap = 1:length(t)
		orb  = 1
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "red" , main = "One orbit" , ticktype = "detailed" )
		
		## Plot4
		snap = 1:length(t)
		orb  = size
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "red" , main = "Another orbit" , ticktype = "detailed" )
	}
}
##}}}

test_rossler = function( plot = FALSE )##{{{
{
	t  = base:::seq( 0 , 50 , length = 1000 )
	size = 1000
	ross = Rossler$new( size = size )
	X   = ross$orbit( t = t )
	
	if( plot )
	{
		
		grDevices::dev.new()
		graphics::par( mfrow = base::c( 2 , 2 ) )
		
		## Plot0
		snap = 1
		orb  = 1:size
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "red" , main = "Initial snapshot" , ticktype = "detailed" )
		
		## Plot1
		snap = length(t)
		orb  = 1:size
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "red" , main = "Final snapshot" , ticktype = "detailed" )
		
		## Plot3
		snap = 1:length(t)
		orb  = 1
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "red" , main = "One orbit" , ticktype = "detailed" )
		
		## Plot4
		snap = 1:length(t)
		orb  = size
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "red" , main = "Another orbit" , ticktype = "detailed" )
	}
}
##}}}

test_lorenz84 = function( plot = FALSE )##{{{
{
	t  = base:::seq( 0 , 2 * 73 , length = 200 * 2 * 73 )
	size = 1000
	l84 = Lorenz84$new( size = size , F = "cyclic" )
	X   = l84$orbit( t = t )
	
	if( plot )
	{
		grDevices::dev.new()
		graphics::par( mfrow = base::c( 2 , 2 ) )
		
		## Plot0
		snap = 400 * 73
		orb  = 1:size
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "orange" , main = "Fall" , ticktype = "detailed" )
		
		## Plot1
		snap = 250 * 73
		orb  = 1:size
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "blue" , main = "Winter" , ticktype = "detailed" )
		
		## Plot2
		snap = 300 * 73
		orb  = 1:size
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "darkgreen" , main = "Spring" , ticktype = "detailed" )
		
		## Plot3
		snap = 350 * 73
		orb  = 1:size
		plot3D::scatter3D( X[snap,orb,1] , X[snap,orb,2] , X[snap,orb,3] , colvar = NULL , col = "red" , main = "Summer" , ticktype = "detailed" )
	}
}
##}}}


test_henon = function( plot = FALSE )##{{{
{
	size = 1000
	henon = Henon$new( size = size )
	t = 1000
	X = henon$orbit(t)
	
	if( plot )
	{
		grDevices::dev.new()
		graphics::par( mfrow = base::c( 2 , 2 ) )
		
		## Plot0
		snap = 1
		orb  = 1:size
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "Initial snapshot" )
		
		## Plot1
		snap = t
		orb  = 1:size
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "Final snapshot" )
		
		## Plot3
		snap = 1:t
		orb  = 1
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "One orbit" )
		
		## Plot4
		snap = 1:t
		orb  = size
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "Another orbit" )
		
	}
}
##}}}

test_ikeda = function( plot = FALSE )##{{{
{
	size = 1000
	ike = Ikeda$new( size = size )
	t = 1000
	X = ike$orbit(t)
	
	if( plot )
	{
		grDevices::dev.new()
		graphics::par( mfrow = base::c( 2 , 2 ) )
		
		## Plot0
		snap = 1
		orb  = 1:size
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "Initial snapshot" )
		
		## Plot1
		snap = t
		orb  = 1:size
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "Final snapshot" )
		
		## Plot3
		snap = 1:t
		orb  = 1
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "One orbit" )
		
		## Plot4
		snap = 1:t
		orb  = size
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "Another orbit" )
		
	}
}
##}}}

test_mira = function( plot = FALSE )##{{{
{
	size = 1000
	mira = Mira$new( size = size )
	t = 1000
	X = mira$orbit(t)
	
	if( plot )
	{
		grDevices::dev.new()
		graphics::par( mfrow = base::c( 2 , 2 ) )
		
		## Plot0
		snap = 1
		orb  = 1:size
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "Initial snapshot" )
		
		## Plot1
		snap = t
		orb  = 1:size
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "Final snapshot" )
		
		## Plot3
		snap = 1:t
		orb  = 1
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "One orbit" )
		
		## Plot4
		snap = 1:t
		orb  = size
		graphics::plot( X[snap,orb,1] , X[snap,orb,2] , col = "red" , main = "Another orbit" )
		
	}
}
##}}}


test_local_dimension = function( plot = FALSE )##{{{
{
	## generate a long orbit
	t  = base:::seq( 0 , 50 , length = 5000 )
	size = 1
	l63 = Lorenz63$new( size = size )
	X0  = l63$orbit( t = base::seq( 0 , 10 , length = 1000 ) )[1000,]
	Y   = l63$orbit( t = t , X0 = X0 )
	
	## Generate a snapshot
	l63 = Lorenz63$new( size = 1000 )
	X   = l63$orbit( t = base::seq( 0 , 20 , length = 1000 ) )[1000,,]
	
	## Find local dimension at "snapshot point" with orbit Y (if Y == NULL, then Y = X)
	ldth = localDimension( X , Y )
	
	if( plot )
	{
		grDevices::dev.new()
		graphics::par( mfrow = base::c( 1 , 2 ) )
		
		## Plot0
		plot3D::scatter3D( X[,1] , X[,2] , X[,3] , colvar = ldth$ld , col = plot3D::ramp.col( base::c("blue","yellow","red") ) , clim = base::c(0.5 , 3.5 ) , main = base::paste( "Local dimension, mean : " , base::round( base::mean(ldth$ld) ,2 ) ) )
		
		## Plot1
		plot3D::scatter3D( X[,1] , X[,2] , X[,3] , colvar = ldth$th , col = plot3D::ramp.col( base::c("blue","yellow","red") ) , main = "Persistence" )
	}
}
##}}}


test_mandelbrot = function( plot = FALSE )##{{{
{
	m = Mandelbrot$new( classic_set = "set0" , nx = 200 , ny = 200 )
	m$run()
	
	if( plot )
	{
		m$plot()
	}
}
##}}}

test_julia = function( plot = FALSE )##{{{
{
	m = Julia$new( classic_set = "set2" , nx = 200 , ny = 200 )
	m$run()
	
	if( plot )
	{
		m$plot()
	}
}
##}}}


run_all_tests = function( plot = FALSE )##{{{
{
	
	test_lorenz63( plot )
	test_rossler( plot )
	test_lorenz84( plot )
	
	test_ikeda( plot )
	test_mira( plot )
	test_henon( plot )
	
	test_local_dimension( plot )
	
	test_mandelbrot( plot )
	test_julia( plot )
}
##}}}



##########
## main ##
##########

#run_all_tests(TRUE)

test_julia(TRUE)


base::cat("Done\n")
