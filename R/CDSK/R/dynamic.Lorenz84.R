
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


