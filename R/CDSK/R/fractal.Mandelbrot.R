
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


###############
## Libraries ##
###############


###############
## Functions ##
###############


## Mandelbrot

#' Mandelbrot
#'
#' Mandelbrot set
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#'        x-coordinates of Mandelbrot set
#' @param y [vector]
#'        y-coordinates of Mandelbrot set
#' @param classic_set [string]
#'        A pre-defined set, values are "set0", "set1" , "set2" , "set3" or "set4".
#' @param nx [integer]
#'        Numbers of step in x-axis
#' @param ny [integer]
#'        Numbers of step in y-axis
#' @param maxit [integer]
#'        Number of iterations to estimate the Mandelbrot set
#' @param new_device [bool]
#'        boolean indicating if grDevices::dev.new() need to be called
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(x,y,classic_set,nx,ny,maxit)}}{This method is used to create object of this class with \code{Mandelbrot}}
#'   \item{\code{run()}}{Estimate the ratio of the Mandelbrot set}
#'   \item{\code{plot()}}{Plot}
#' }
#' @examples
#' m = Mandelbrot$new( classic_set = "set0" , nx = 200 , ny = 200 )
#' m$run()
#'
#' @export
Mandelbrot = R6::R6Class( "Mandelbrot" ,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	ratio = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( x = NULL , y = NULL , classic_set = NULL , nx = NULL , ny = NULL , maxit = 200 )
	{
		private$x = x
		private$y = y
		private$nx = nx
		private$ny = ny
		private$maxit = maxit
		
		if( !is.null(classic_set) && ( !is.null(nx) || !is.null(ny) ) )
		{
			private$nx = if( is.null(nx) ) ny else nx
			private$ny = if( is.null(ny) ) nx else ny
			
			if(      classic_set == "set0" ) private$set_set0()
			else if( classic_set == "set1" ) private$set_set1()
			else if( classic_set == "set2" ) private$set_set2()
			else if( classic_set == "set3" ) private$set_set3()
			else if( classic_set == "set4" ) private$set_set4()
		}
		else if( !is.null(x) && !is.null(y) )
		{
			private$nx = length(x)
			private$ny = length(y)
		}
		else
		{
			print("Error, Need (x,y) set or (classic_set,nx [,ny]) set!" )
		}
	},
	
	
	#############
	## Methods ##
	#############
	
	run = function()
	{
		grid = as.matrix( base::expand.grid( private$y , private$x ) )
		c = complex( 1 , grid[,2] , grid[,1] )
		z = numeric( private$nx * private$ny )
		ratio = numeric( private$nx * private$ny )
		
		for( i in 1:private$maxit )
		{
			z = z^2 + c
			idx = base::abs(z) < 2
			ratio[idx] = ratio[idx] + 1
		}
		
		self$ratio = matrix( ratio , nrow = private$nx , ncol = private$ny ) / private$maxit
	},
	
	plot = function( new_device = FALSE )
	{
		if(new_device)
			grDevices::dev.new()
		graphics::image( self$ratio , col = base::rev(heat.colors(12)) )
	}
	
	),
	
	
	private = list(
	
	###############
	## arguments ##
	###############
	
	x = NULL,
	y = NULL,
	nx    = NULL,
	ny    = NULL,
	maxit = NULL,
	
	
	#############
	## Methods ##
	#############
	
	set_set0 = function()
	{
		private$x = base::seq( -0.62 , -0.42 , length = private$nx )
		private$y = base::seq( -0.7  , -0.5  , length = private$ny )
	},
	
	set_set1 = function()
	{
		private$x = base::seq( -0.6 , -0.595 , length = private$nx )
		private$y = base::seq( -0.665  , -0.66  , length = private$ny )
	},
	
	set_set2 = function()
	{
		private$x = base::seq( -0.57 , -0.56 , length = private$nx )
		private$y = base::seq( -0.646  , -0.636  , length = private$ny )
	},
	
	set_set3 = function()
	{
		private$x = base::seq( -0.5663 , -0.5656 , length = private$nx )
		private$y = base::seq( -0.6394  , -0.6387  , length = private$ny )
	},
	
	set_set4 = function()
	{
		private$x = base::seq( -1.5 , -0.5 , length = private$nx )
		private$y = base::seq( -0.5  , 0.5  , length = private$ny )
	}
	
	)
)

