
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


## DiscDynSystem {{{

#' DiscDynSystem
#'
#' @description
#' Base class to define Differential Dynamical system.
#'
#' @details
#' Do not use it!!!
#'
#' @examples
#' ## No example because you should not use this class!
#'
#' @export
DiscDynSystem = R6::R6Class( "DiscDynSystem" ,
	
	inherit = CDSK::DynamicalSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	
	#################
	## Constructor ##
	#################
	
	#' @description
    #' Create a new DiscDynSystem object.
	#' @param dim    [integer] Dimension of the phase space of the dynamical
	#'        system
	#' @param size   [integer] Number of initial condition simultaneously solved
	#' @param bounds [matrix]  Bounds of phase space where initial condition can
	#'        be drawn
	#' @return A new `DiscDynSystem` object.
	initialize = function( dim , size , bounds )
	{
		super$initialize( dim , size , bounds )
	}
	
	),
	
	private = list(
	
	#############
	## Methods ##
	#############
	
	
	solver = function( t , X0 )
	{
		n = if( length(t) == 1 ) as.integer(t) else length(t)
		X = base::matrix( NA , n + 1 , self$size * self$dim )
		X[1,] = X0
		for( i in 1:n )
		{
			X[i+1,] = private$equation( NULL , X[i,] , NULL )
		}
		return( X )
	}
	
	
	)
)
##}}}
