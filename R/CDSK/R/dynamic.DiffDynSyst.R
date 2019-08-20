
###############
## Libraries ##
###############


###############
## Functions ##
###############


## DiffDynSystem {{{

#' DiffDynSystem
#'
#' Base class to define Differential Dynamical system, do not use it!!!
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom deSolve ode
#'
#' @param dim [integer]
#'        Dimension of the phase space of the dynamical system
#' @param size [integer]
#'        Number of initial condition simultaneously solved
#' @param bounds  [matrix]
#'        Bounds of phase space where initial condition can be drawn
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
#'   \item{\code{new(dim,size,bounds)}}{This method is used to create object of this class with \code{DiffDynSystem}}
#' }
#' @examples
#' ## No example because you should not use this class!
#'
#' @export
DiffDynSystem = R6::R6Class( "DiffDynSystem" ,
	
	inherit = CDSK::DynamicalSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	
	#################
	## Constructor ##
	#################
	
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
		X = deSolve::ode( X0 , t , private$equation , NULL , method = "rk4" )[,-1]
		return( X )
	}
	
	
	)
)
##}}}


