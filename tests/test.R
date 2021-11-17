
## Copyright(c) 2021 Yoann Robin
## 
## This file is part of CDSK.
## 
## CDSK is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## CDSK is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with CDSK.  If not, see <https://www.gnu.org/licenses/>.


base::rm( list = base::ls() )

###############
## Libraries ##
###############

library(devtools)
load_all( "../R/CDSK" )
#library(CDSK)
library(plot3D)


################
## Plot tools ##
################

PlotTools = R6::R6Class( "PlotTools" , ##{{{
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	os = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
		self$os = self$get_os()
	},
	
	
	#############
	## Methods ##
	#############
	
	get_os = function()
	{
		sysinf = base::Sys.info()
		if( !is.null(sysinf) )
		{
			os = sysinf['sysname']
			if( os == 'Darwin' ) os = "osx"
		}
		else
		{
			## mystery machine
			os = .Platform$OS.type
			if( base::grepl( "^darwin"   , R.version$os ) ) os = "osx"
			if( base::grepl( "linux-gnu" , R.version$os ) ) os = "linux"
		}
		invisible(tolower(os))
	},
	
	new_screen = function()
	{
		if( self$os == "osx" )
		{
			grDevices::quartz()
		}
		if( self$os == "linux" )
		{
			grDevices::X11()
		}
	},
	
	wait = function()
	{
		while( base::names(grDevices::dev.cur()) !='null device' ) base::Sys.sleep(1)
	}
	
	)
)
##}}}

plt = PlotTools$new()


###############
## Functions ##
###############

abind = function( X0 , X1 )##{{{
{
	X = array( NA , dim = base::c(dim(X0),2) )
	X[,,1] = X0
	X[,,2] = X1
	return(X)
}
##}}}


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
		plt$new_screen()
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
		
		plt$new_screen()
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
		plt$new_screen()
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
		plt$new_screen()
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
		plt$new_screen()
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
		plt$new_screen()
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


test_dynamical_local_indexes = function()##{{{
{
	## Data
	l63 = CDSK::Lorenz63$new( size = 1000 )
	X01 = l63$orbit( base::seq( 0 , 50 , length = 1000 ) )[1000,,]
	l63 = CDSK::Lorenz63$new()
	y0  = l63$orbit( base::seq( 0 , 50 , length = 1000 ) )[1000,]
	y1  = y0 + stats::rnorm( mean = 0 , sd = 0.1 , n = 3 )
	Y0  = l63$orbit( base::seq( 0 , 50 , length = 1000 ) , y0 )
	Y1  = l63$orbit( base::seq( 0 , 50 , length = 1000 ) , y1 )
	
	X = abind( X01 , X01 )
	Y = abind(  Y0 ,  Y1 )
	
	## Args
	dli = dynamical_local_indexes( X , Y , return_shape = TRUE )

}
##}}}


test_mandelbrot = function( plot = FALSE )##{{{
{
	m = Mandelbrot$new( classic_set = "set0" , nx = 200 , ny = 200 )
	m$run()
	
	if( plot )
	{
		plt$new_screen()
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
		plt$new_screen()
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
	
	test_dynamical_local_indexes()
	
	test_mandelbrot( plot )
	test_julia( plot )
}
##}}}


##########
## main ##
##########

## Read command line arguments and run (or not) tests
##================================================{{{

args = commandArgs( trailingOnly = TRUE )
args_verbose = FALSE
args_run     = FALSE
if( length(args) > 0 )
{
	for( a in args )
	{
		if( a == "-r" || a == "--run-all-tests" )
			args_run = TRUE
		if( a == "-v" || a == "--verbose" )
			args_verbose = TRUE
	}
}

if( args_run )
	run_all_tests(args_verbose)

##}}}

plt$wait()

base::cat("Done\n")
