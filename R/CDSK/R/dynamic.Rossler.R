
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


## Rossler {{{

#' Rossler
#'
#' @description
#' Rossler dynamical system
#'
#' @details
#' Continuous dynamical system.
#'
#' @references O. E. Rossler, « An equation for continuous chaos », Phys. Rev.
#'             Lett. A, vol. 57, np5, 1976, p. 397-398
#'
#' @examples
#' ross = CDSK::Rossler$new( size = 200 )
#' t = base::seq( 0 , 100 , 0.005 )
#' X = ross$orbit(t) 
#' ## X is an array with dim = (length(t),200,3)
#' ## Each X[,i,] is an orbit
#' ## Each X[i,,] is a snapshot
#'
#' @export
Rossler = R6::R6Class( "Rossler" ,
	
	inherit = CDSK::DiffDynSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	#' @field a [float] Default = 0.1
	a = 0.1,
	#' @field b [float] Default = 0.1
	b = 0.1,
	#' @field c [float] Default = 14
	c = 14.,
	
	#################
	## Constructor ##
	#################
	
	#' @description
    #' Create a new Rossler object.
	#' @param a [float] Default = 0.1
	#' @param b [float] Default = 0.1
	#' @param c [float] Default = 14
	#' @param size [integer] Number of initial condition simultaneously solved
	#' @return A new `Rossler` object.
	initialize = function( a = 0.1 , b = 0.1 , c = 14 , size = 1 )
	{
		super$initialize( 3 , size , base::matrix( base::c( -20 , -20 , 0 , 20 , 20 , 35 ) , nrow = 3 , ncol = 2 ) )
		self$a = a
		self$b = b
		self$c = c
	}
	
	
	#############
	## Methods ##
	#############
	
	),
	
	private = list(
	
	###############
	## arguments ##
	###############
	
	#############
	## Methods ##
	#############
	
	equation = function( t , X , par = NULL )
	{
		dX = numeric( length(X) )
		dX[private$i[[1]]] = - X[private$i[[2]]] - X[private$i[[3]]]
		dX[private$i[[2]]] = X[private$i[[1]]] + self$a * X[private$i[[2]]]
		dX[private$i[[3]]] = self$b + X[private$i[[3]]] * ( X[private$i[[1]]] - self$c )
		
		return(list(dX))
	}
	
	)
)
##}}}

