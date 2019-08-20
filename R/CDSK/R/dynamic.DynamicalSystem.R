
###############
## Libraries ##
###############


###############
## Functions ##
###############

## DynamicalSystem {{{

#' DynamicalSystem
#'
#' Base class to define Dynamical system, do not use it!!!
#'
#' @docType class
#' @importFrom R6 R6Class
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
#'   \item{\code{new(dim,size,bounds)}}{This method is used to create object of this class with \code{DynamicalSystem}}
#'   \item{\code{orbit(t,X0)}}{Compute the orbit along t starting at X0. If X0 is NULL, it is randomly drawn by randomIC()}
#'   \item{\code{randomIC()}}{Return a random initial condition}
#' }
#' @examples
#' ## No example because you should not use this class!
#'
#' @export
DynamicalSystem = R6::R6Class( "DynamicalSystem" ,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	dim    = 0,
	size   = 0,
	bounds = NULL,
	
	
	#################
	## Constructor ##
	#################
	
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
	
	randomIC = function()
	{
		return( stats::runif( self$dim * self$size , min = self$bounds[,1] , max = self$bounds[,2] ) )
#		return( base::apply( self$bounds , 1 , function(X) { return( stats::runif( n = self$size , min = X[1] , max = X[2] ) ) } ) )
	},
	
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


