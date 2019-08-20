

###############
## Libraries ##
###############


###############
## Functions ##
###############

## Mira {{{

#' Mira
#'
#' Mira dynamical system
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param a [float]
#'		  Default = -0.48
#' @param b [float]
#'        Default = 0.93
#' @param size [integer]
#'        Number of initial condition simultaneously solved
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
#'   \item{\code{new(a,b,size)}}{This method is used to create object of this class with \code{Mira}}
#'   \item{\code{orbit(t,X0)}}{Compute the orbit along t starting at X0. If X0 is NULL, it is randomly drawn by randomIC()}
#'   \item{\code{randomIC()}}{Return a random initial condition}
#' }
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
	
	a = -0.48,
	b = 0.93,
	
	
	#################
	## Constructor ##
	#################
	
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

