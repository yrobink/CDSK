
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

## DynamicalSystem {{{

#' DynamicalSystem
#'
#' @description
#' Base class to define Dynamical system
#'
#' @details
#' Do not use it!!!
#'
#' @examples
#' ## No example because you should not use this class!
#'
#' @export
DynamicalSystem = R6::R6Class( "DynamicalSystem" ,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	#' @field dim    [integer] Dimension of the phase space of the dynamical
	#'        system
	dim    = 0,
	#' @field size   [integer] Number of initial condition simultaneously solved
	size   = 0,
	#' @field bounds [matrix]  Bounds of phase space where initial condition can
	#'        be drawn
	bounds = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	#' @description
    #' Create a new DynamicalSystem object.
	#' @param dim    [integer] Dimension of the phase space of the dynamical
	#'        system
	#' @param size   [integer] Number of initial condition simultaneously solved
	#' @param bounds [matrix]  Bounds of phase space where initial condition can
	#'        be drawn
	#' @return A new `DynamicalSystem` object.
	initialize = function( dim , size = 1 , bounds = NULL )
	{
		self$dim    = dim
		self$size   = size
		self$bounds = bounds
		
		private$i = list()
		for( i in 1:self$dim )
		{
			private$i[[i]] = seq( i , self$dim * self$size , self$dim )
		}
	},
	
	
	#############
	## Methods ##
	#############
	
	#' @description
    #' Draw a random initial condition.
    #' @return [vector] A random initial condition.
	randomIC = function()
	{
		return( stats::runif( self$dim * self$size , min = self$bounds[,1] , max = self$bounds[,2] ) )
	},
	
	#' @description
    #' Compute an orbit of the dynamical system.
	#' @param t  [vector] Time to integrate the dynamical system
	#' @param X0 [vector or NULL] Vector of initial condition of size size*dim,
	#'        if NULL a random IC is drawn in bounds
    #' @return   [matrix] The orbit.
	orbit = function( t , X0 = NULL )
	{
		X0 = if( is.null(X0) ) self$randomIC() else as.vector(base::t(X0))
		X  = private$solver( t , X0 )
		if( self$size > 1 )
		{
			X  = base::aperm( base::array( X , base::c( dim(X)[1] , self$dim , self$size ) ) , perm = base::c( 1 , 3 , 2 ) )
		}
		return( X )
	}
	
	
	),
	
	private = list(
	
	###############
	## arguments ##
	###############
	
	i = NULL,
	
	
	#############
	## Methods ##
	#############
	
	solver = function( t , X0 )
	{},
	
	equation = function( t , X , par = NULL )
	{}
	
	)
)
##}}}


