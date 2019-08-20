
###############
## Libraries ##
###############


###############
## Functions ##
###############


## Rossler {{{

#' Rossler
#'
#' Rossler dynamical system
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param a [float]
#'        Default = 0.1
#' @param b [float]
#'        Default = 0.1
#' @param c [float]
#'        Default = 14
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
#'   \item{\code{new(a,b,c,size)}}{This method is used to create object of this class with \code{Rossler}}
#'   \item{\code{orbit(t,X0)}}{Compute the orbit along t starting at X0. If X0 is NULL, it is randomly drawn by randomIC()}
#'   \item{\code{randomIC()}}{Return a random initial condition}
#' }
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
	
	a = 0.1,
	b = 0.1,
	c = 14.,
	
	#################
	## Constructor ##
	#################
	
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

