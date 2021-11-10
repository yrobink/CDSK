

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the CDSK (Chaotic        ##
## Dynamical System Kit) library. This library makes it possible                ##
## to generate classic (continuous and discrete) attractors, generate the       ##
## Mandelbrot and Julia set, and fit the local dimension.                       ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## CDSK (Chaotic Dynamical System Kit). Cette librairie permet de générer les   ##
## attracteurs classiques (discret comme continue), de générer l'ensemble de    ##
## Julia et de Mandelbrot et d'estimer les dimensions locales.                  ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################


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
		plt$new_screen()
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
	
	test_local_dimension( plot )
	
	test_mandelbrot( plot )
	test_julia( plot )
}
##}}}

test_dynamical_local_indexes = function()
{
	## Data
	l63 = CDSK::Lorenz63$new( size = 1000 )
	X01 = l63$orbit( base::seq( 0 , 100 , length = 10000 ) )[10000,,]
	l63 = CDSK::Lorenz63$new()
	y0  = l63$orbit( base::seq( 0 , 100 , length = 10000 ) )[10000,]
	y1  = y0 + stats::rnorm( mean = 0 , sd = 0.1 , n = 3 )
	Y0  = l63$orbit( base::seq( 0 , 100 , length = 10000 ) , y0 )
	Y1  = l63$orbit( base::seq( 0 , 100 , length = 10000 ) , y1 )
	
	X = abind( X01 , X01 )
	Y = abind(  Y0 ,  Y1 )
	
	## Args
	dli = dynamical_local_indexes( X , Y )

}


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
