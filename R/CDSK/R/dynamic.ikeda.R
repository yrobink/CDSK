
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

## Ikeda {{{

#' Ikeda
#'
#' @description
#' Ikeda dynamical system
#'
#' @details
#' Discrete dynamical system.
#'
#' @examples
#' ike = CDSK::Ikeda$new( size = 200 )
#' X = ike$orbit(1000) 
#' ## X is an array with shape = (1000,200,2)
#' ## Each X[,i,] is an orbit
#' ## Each X[i,,] is a snapshot
#'
#' @export
Ikeda = R6::R6Class( "Ikeda" , 
	
	inherit = CDSK::DiscDynSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	#' @field R  [float] Default = 1
	R = 1.,
	#' @field C1 [float] Default = 0.4
	C1 = 0.4,
	#' @field C2 [float] Default = 0.9
	C2 = 0.9,
	#' @field C3 [float] Default = 6.
	C3 = 6.,
	
	#################
	## Constructor ##
	#################
	
	#' @description
    #' Create a new Ikeda object.
	#' @param R  [float] Default = 1
	#' @param C1 [float] Default = 0.4
	#' @param C2 [float] Default = 0.9
	#' @param C3 [float] Default = 6.
	#' @param size [integer] Number of initial condition simultaneously solved
	#' @return A new `Ikeda` object.
	initialize = function( R = 1. , C1 = 0.4 , C2 = 0.9 , C3 = 6. , size = 1 )
	{
		super$initialize( 2 , size , base::matrix( base::c( 0 , -1 , 1 , 0 ) , nrow = 2 , ncol = 2 ) )
		self$R  = R
		self$C1 = C1
		self$C2 = C2
		self$C3 = C3
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
		tau = self$C1 - self$C3 / ( 1. + X[private$i[[1]]]^2 + X[private$i[[2]]]^2 )
		Xnext[private$i[[1]]] = self$R + self$C2 * ( X[private$i[[1]]] * base::cos(tau) - X[private$i[[2]]] * base::sin(tau) )
		Xnext[private$i[[2]]] = self$C2 * ( X[private$i[[2]]] * base::cos(tau) + X[private$i[[1]]] * base::sin(tau) )
		return(Xnext)
	}
	
	)
)
##}}}

