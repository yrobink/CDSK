
################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2019                                                ##
##                                                                            ##
## yoann.robin.k@gmail.com                                                    ##
##                                                                            ##
## This software is a computer program that is part of the CDSK (Chaotic      ##
## Dynamical System Kit) library. This library makes it possible              ##
## to generate classic (continuous and discrete) attractors, generate the     ##
## Mandelbrot and Julia set, and fit the local dimension.                     ##
##                                                                            ##
## This software is governed by the CeCILL-C license under French law and     ##
## abiding by the rules of distribution of free software.  You can  use,      ##
## modify and/ or redistribute the software under the terms of the CeCILL-C   ##
## license as circulated by CEA, CNRS and INRIA at the following URL          ##
## "http://www.cecill.info".                                                  ##
##                                                                            ##
## As a counterpart to the access to the source code and  rights to copy,     ##
## modify and redistribute granted by the license, users are provided only    ##
## with a limited warranty  and the software's author,  the holder of the     ##
## economic rights,  and the successive licensors  have only  limited         ##
## liability.                                                                 ##
##                                                                            ##
## In this respect, the user's attention is drawn to the risks associated     ##
## with loading,  using,  modifying and/or developing or reproducing the      ##
## software by the user in light of its specific status of free software,     ##
## that may mean  that it is complicated to manipulate,  and  that  also      ##
## therefore means  that it is reserved for developers  and  experienced      ##
## professionals having in-depth computer knowledge. Users are therefore      ##
## encouraged to load and test the software's suitability as regards their    ##
## requirements in conditions enabling the security of their systems and/or   ##
## data to be ensured and,  more generally, to use and operate it in the      ##
## same conditions as regards security.                                       ##
##                                                                            ##
## The fact that you are presently reading this means that you have had       ##
## knowledge of the CeCILL-C license and that you accept its terms.           ##
##                                                                            ##
################################################################################
################################################################################

################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2019                                                ##
##                                                                            ##
## yoann.robin.k@gmail.com                                                    ##
##                                                                            ##
## Ce logiciel est un programme informatique faisant partie de la librairie   ##
## CDSK (Chaotic Dynamical System Kit). Cette librairie permet de générer les ##
## attracteurs classiques (discret comme continue), de générer l'ensemble de  ##
## Julia et de Mandelbrot et d'estimer les dimensions locales.                ##
##                                                                            ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et  ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez    ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions     ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA   ##
## sur le site "http://www.cecill.info".                                      ##
##                                                                            ##
## En contrepartie de l'accessibilité au code source et des droits de copie,  ##
## de modification et de redistribution accordés par cette licence, il n'est  ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le     ##
## titulaire des droits patrimoniaux et les concédants successifs.            ##
##                                                                            ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques      ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au      ##
## développement et à la reproduction du logiciel par l'utilisateur étant     ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à      ##
## manipuler et qui le réserve donc à des développeurs et des professionnels  ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les    ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du     ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la       ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,   ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         ##
##                                                                            ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez     ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les  ##
## termes.                                                                    ##
##                                                                            ##
################################################################################
################################################################################


###############
## Libraries ##
###############


###############
## Functions ##
###############

## Lorenz84TimeForcing {{{

#' Lorenz84TimeForcing
#'
#' @description
#' Lorenz84 time forcing.
#'
#' @details
#' Can be constant, cyclic (period = 73), linear, or cyclic and linear.
#'
#' @param t [vector]
#'        Time to evaluate the forcing
#'
#' @examples
#' ## No example, used by Lorenz84 model
#'
#' @export
Lorenz84TimeForcing = R6::R6Class( "Lorenz84TimeForcing" ,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	#' @field tcc [float] Time of Climate Change, default = 100 * 73 (100 years)
	tcc = 0,
	
	
	#################
	## Constructor ##
	#################
	
	#' @description
    #' Create a new Lorenz84TimeForcing object.
	#' @param tcc [float] Time of Climate Change, default = 100 * 73 (100 years)
	#' @return A new `Lorenz84TimeForcing` object.
	initialize = function( tcc = 100 * 73 )
	{
		self$tcc = tcc
	},
	
	
	#############
	## Methods ##
	#############
	
	#' @description
    #' Constant forcing.
	#' @param t [vector] Time to evaluate the forcing
	#' @return [vector] value of forcing.
	constant = function(t)
	{
		return(6)
	},
	
	#' @description
    #' Cyclic forcing, similar to seasonnal cycle
	#' @param t [vector] Time to evaluate the forcing
	#' @return [vector] value of forcing.
	cyclic = function(t)
	{
		return( 9.5 + 2 * base::sin( 2 * base::pi / 73 * t ) )
	},
	
	#' @description
    #' Linear forcing, similar to climate change
	#' @param t [vector] Time to evaluate the forcing
	#' @return [vector] value of forcing.
	linear = function( t )
	{
		if( t < self$tcc )
		{
			return(0)
		}
		else
		{
			return( - 2 * ( t - self$tcc ) / self$tcc )
		}
	}
	
	)
)
##}}}

## Lorenz84 {{{

#' Lorenz84
#'
#' @description
#' Lorenz84 dynamical system
#'
#' @details
#' Continuous dynamical system.
#'
#' @references E. N. Lorenz, "Irregularity : a Fundamental property of the
#'             atmosphere", In Tellus A, vol. 36, n. 2, p. 98-110 (1984)
#'
#' @references E. N. Lorenz, "Can chaos and intransitivity lead to interannual
#'             variability?", In Tellus A, vol.42, n. 3, p. 378-389 (1990)
#'
#' @references G. Drotos and al, “Probabilistic concepts in a changing climate :
#'             a snapshot attractor picture”. In : Jour. Clim., vol. 28, n. 8,
#'             p. 3275–3288 (2015)
#'
#' @examples
#' l84 = CDSK::Lorenz84$new( size = 200 , F = "cyclic" )
#' t = base::seq( 0 , 100 , 0.005 )
#' X = l84$orbit(t) 
#' ## X is an array with dim = (length(t),200,3)
#' ## Each X[,i,] is an orbit
#' ## Each X[i,,] is a snapshot
#'
#' @export
Lorenz84 = R6::R6Class( "Lorenz84" , 
	
	inherit = CDSK::DiffDynSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	#' @field a [float] Default = 0.25
	a = 0.25,
	#' @field b [float] Default = 4.
	b = 4.,
	#' @field G [float] Default = 1.
	G = 1.,
	#' @field F [callable or string] Time forcing.
	F = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	#' @description
    #' Create a new Lorenz84 object.
	#' @param a [float] Default = 0.25
	#' @param b [float] Default = 4.
	#' @param G [float] Default = 1.
	#' @param F [callable or string] Time forcing, if string:
	#'        => "constant" use Lorenz84TimeForcing$constant
	#'        => "cyclic" use Lorenz84TimeForcing$cyclic
	#'        => "linear" use Lorenz84TimeForcing$linear
	#'        => "cyclic-linear" use cyclic + linear
	#' @param size [integer] Number of initial condition simultaneously solved
	#' @return A new `Lorenz84` object.
	initialize = function( a = 0.25 , b = 4. , G = 1. , F = NULL , size = 1 )
	{
		super$initialize( 3 , size , base::matrix( base::c( -1 , -3 , -3 , 3 , 3 , 3 ) , nrow = 3 , ncol = 2 ) )
		self$a = a
		self$b = b
		self$G = G
		if( is.function(F) )
		{
			self$F = F
		}
		else if( is.character(F) )
		{
			private$forcing = CDSK::Lorenz84TimeForcing$new()
			if( F == "cyclic" )
			{
				self$F = private$forcing$cyclic
			}
			else if( F == "linear" )
			{
				self$F = private$forcing$linear
			}
			else if( F == "cyclic-linear" )
			{
				self$F = function(t) { return( private$forcing$cyclic(t) + private$forcing$linear(t) ) }
			}
		}
		else
		{
			private$forcing = CDSK::Lorenz84TimeForcing$new()
			self$F = private$forcing$constant
		}
	}
	
	
	#############
	## Methods ##
	#############
	
	),
	
	private = list(
	
	###############
	## arguments ##
	###############
	
	forcing = NULL,
	
	
	#############
	## Methods ##
	#############
	
	equation = function( t , X , par = NULL )
	{
		dX = numeric( length(X) )
		dX[private$i[[1]]] = - X[private$i[[2]]]^2 - X[private$i[[3]]]**2 - self$a * X[private$i[[1]]] + self$a * self$F(t)
		dX[private$i[[2]]] = X[private$i[[1]]] * X[private$i[[2]]] - self$b * X[private$i[[1]]] * X[private$i[[3]]] - X[private$i[[2]]] + self$G
		dX[private$i[[3]]] = X[private$i[[1]]] * X[private$i[[3]]] + self$b * X[private$i[[1]]] * X[private$i[[2]]] - X[private$i[[3]]]
		return(list(dX))
	}
	
	)
)
##}}}


