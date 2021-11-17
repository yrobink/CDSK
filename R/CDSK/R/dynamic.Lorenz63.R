
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

## Lorenz63 {{{

#' Lorenz63
#'
#' @description
#' Lorenz (1963) dynamical system
#'
#' @details
#' Continuous dynamical system.
#'
#' @references E. N. Lorenz, « Deterministic nonperiodic flow », J. Atmos. Sci.,
#'             vol. 20, no 2, 1963, p. 130-141
#'
#' @examples
#' l63 = CDSK::Lorenz63$new( size = 200 )
#' t = base::seq( 0 , 100 , 0.005 )
#' X = l63$orbit(t) 
#' ## X is an array with dim = (length(t),200,3)
#' ## Each X[,i,] is an orbit
#' ## Each X[i,,] is a snapshot
#'
#' @export
Lorenz63 = R6::R6Class( "Lorenz63" , 
	
	inherit = CDSK::DiffDynSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	#' @field s [float] Number of Prandtl, default = 10
	s = 10,
	#' @field r [float] Number of Rayleigh, default = 28
	r = 28,
	#' @field b [float] Ratio of critical value, default = 2.667
	b = 2.667,
	
	#################
	## Constructor ##
	#################
	
	#' @description
    #' Create a new Lorenz63 object.
	#' @param s [float] Number of Prandtl, default = 10
	#' @param r [float] Number of Rayleigh, default = 28
	#' @param b [float] Ratio of critical value, default = 2.667
	#' @param size [integer] Number of initial condition simultaneously solved
	#' @return A new `Lorenz63` object.
	initialize = function( s = 10 , r = 28 , b = 2.667 , size = 1 )
	{
		super$initialize( 3 , size , base::matrix( base::c( -20 , -20 , 0 , 20 , 20 , 40 ) , nrow = 3 , ncol = 2 ) )
		self$s = s
		self$r = r
		self$b = b
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
		dX[private$i[[1]]] = self$s * ( X[private$i[[2]]] - X[private$i[[1]]] )
		dX[private$i[[2]]] = self$r * X[private$i[[1]]] - X[private$i[[2]]] - X[private$i[[1]]] * X[private$i[[3]]]
		dX[private$i[[3]]] = X[private$i[[1]]] * X[private$i[[2]]] - self$b * X[private$i[[3]]]
		
		return(list(dX))
	}
	
	)
)
##}}}

