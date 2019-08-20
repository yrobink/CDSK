
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

## Lorenz84TimeForcing {{{

#' Lorenz84TimeForcing
#'
#' Lorenz84 time forcing, can be constant, cyclic (period = 73), linear, or cyclic and linear.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param tcc [float]
#'        Time of Climate Change, default = 100 * 73 (100 years)
#' @param t [vector]
#'        Time to evaluate the forcing
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(tcc)}}{This method is used to create object of this class with \code{Lorenz84TimeForcing}}
#'   \item{\code{constant(t)}}{Constant time forcing, fixed at 6}
#'   \item{\code{cyclic(t)}}{Cyclic time forcing, with period fixed at 73, varying between 7.5 and 11.5}
#'   \item{\code{linear(t)}}{Linear forcing, 0 before tcc, and decreasing linearly after tcc.}
#' }
#' @examples
#' ## No example, used by Lorenz84 model
#'
#' @export
Lorenz84TimeForcing = R6::R6Class( "Lorenz84TimeForcing" ,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	tcc = 0,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( tcc = 100 * 73 )
	{
		self$tcc = tcc
	},
	
	
	#############
	## Methods ##
	#############
	
	constant = function(t)
	{
		return(6)
	},
	
	cyclic = function(t)
	{
		return( 9.5 + 2 * base::sin( 2 * base::pi / 73 * t ) )
	},
	
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
#' Lorenz84 dynamical system
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param a [float]
#'        Default = 0.25
#' @param b [float]
#'        Default = 4.
#' @param G [float]
#'        Default = 1.
#' @param F [callable or string]
#'        Time forcing, if string:
#'        => "constant" use Lorenz84TimeForcing$constant
#'        => "cyclic" use Lorenz84TimeForcing$cyclic
#'        => "linear" use Lorenz84TimeForcing$linear
#'        => "cyclic-linear" use Lorenz84TimeForcing$cyclic + Lorenz84TimeForcing$linear
#' @param size [integer]
#'        Number of initial condition simultaneously solved
#' @param t [vector]
#'        Time to integrate the dynamical system
#' @param X0 [vector or NULL]
#'        Vector of initial condition of size size*dim, if NULL a random IC is drawn in bounds
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(a,b,G,F,size)}}{This method is used to create object of this class with \code{Lorenz84}}
#'   \item{\code{orbit(t,X0)}}{Compute the orbit along t starting at X0. If X0 is NULL, it is randomly drawn by randomIC()}
#'   \item{\code{randomIC()}}{Return a random initial condition}
#' }
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
	
	a = 0.25,
	b = 4.,
	G = 1.,
	F = NULL,
	
	
	#################
	## Constructor ##
	#################
	
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


