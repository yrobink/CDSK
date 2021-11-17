
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

## Henon {{{

#' Henon
#'
#' @description
#' Henon dynamical system
#'
#' @details
#' Discrete dynamical system.
#'
#' @examples
#' henon = CDSK::Henon$new( size = 200 )
#' X = henon$orbit(1000) 
#' ## X is an array with shape = (1000,200,2)
#' ## Each X[,i,] is an orbit
#' ## Each X[i,,] is a snapshot
#'
#' @export
Henon = R6::R6Class( "Henon" , 
	
	inherit = CDSK::DiscDynSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	#' @field a [float] Parameter of Henon system
	a = 1.4,
	#' @field b [float] Parameter of Henon system
	b = 0.3,
	
	
	#################
	## Constructor ##
	#################
	
	#' @description
    #' Create a new Henon object.
	#' @param a    [float]   Default = 1.4
	#' @param b    [float]   Default = 0.3
	#' @param size [integer] Number of initial condition simultaneously solved
	#' @return A new `Henon` object.
	initialize = function( a = 1.4 , b = 0.3 , size = 1 )
	{
		super$initialize( 2 , size , base::matrix( base::c( 0 , 0 , 0.5 , 0.5 ) , nrow = 2 , ncol = 2 ) )
		self$a = a
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
		Xnext = numeric(length(X))
		Xnext[private$i[[1]]] = X[private$i[[2]]] + 1 - self$a * X[private$i[[1]]]^2
		Xnext[private$i[[2]]] = self$b * X[private$i[[1]]]
		return(Xnext)
	}
	
	)
)
##}}}

