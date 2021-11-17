
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

## Mira {{{

#' Mira
#'
#' @description
#' Mira dynamical system
#'
#' @details
#' Discrete dynamical system.
#'
#' @examples
#' mira = CDSK::Mira$new( size = 200 )
#' X = mira$orbit(1000) 
#' ## X is an array with shape = (1000,200,2)
#' ## Each X[,i,] is an orbit
#' ## Each X[i,,] is a snapshot
#'
#' @export
Mira = R6::R6Class( "Mira" , 
	
	inherit = CDSK::DiscDynSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	#' @field a [float] Default = -0.48
	a = -0.48,
	#' @field b [float] Default = 0.93
	b = 0.93,
	
	
	#################
	## Constructor ##
	#################
	
	#' @description
    #' Create a new Mira object.
	#' @param a [float] Default = -0.48
	#' @param b [float] Default = 0.93
	#' @param size [integer] Number of initial condition simultaneously solved
	#' @return A new `Mira` object.
	initialize = function( a = -0.48 , b = 0.93 , size = 1 )
	{
		super$initialize( 2 , size , base::matrix( base::c( 3.8 , -0.2 , 4.2 , 0.2 ) , nrow = 2 , ncol = 2 ) )
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
		Xnext[private$i[[1]]] = self$b * X[private$i[[2]]] + private$F( X[private$i[[1]]] )
		Xnext[private$i[[2]]] = - X[private$i[[1]]] + private$F( Xnext[private$i[[1]]] )
		
		return(Xnext)
	},
	
	F = function( x )
	{
		return( self$a * x + 2 * ( 1 - self$a ) * x^2 / ( 1 + x^2 ) )
	}
	
	)
)
##}}}

