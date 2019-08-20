

###############
## Libraries ##
###############


###############
## Functions ##
###############

## Lorenz63 {{{

#' Lorenz63
#'
#' Lorenz (1963) dynamical system
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom deSolve ode
#'
#' @param s [float]
#'        Number of Prandtl, default = 10
#' @param r [float]
#'        Number of Rayleigh, default = 28
#' @param b [float]
#'        Ratio of critical value, default = 2.667
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
#'   \item{\code{new(s,r,b,size)}}{This method is used to create object of this class with \code{Lorenz63}}
#'   \item{\code{orbit(t,X0)}}{Compute the orbit along t starting at X0. If X0 is NULL, it is randomly drawn by randomIC()}
#'   \item{\code{randomIC()}}{Return a random initial condition}
#' }
#' @examples
#' l63 = CDSK::Lorenz63$new( size = 200 )
#' t = base::seq( 0 , 100 , 0.005 )
#' X = l63$orbit(t) 
#' ## X is an array with dim = (length(t),3,200)
#' ## Each X[i,,] is an orbit
#' ## Each X[,,i] is a snapshot
#'
#' @export
Lorenz63 = R6::R6Class( "Lorenz63" , 
	
	inherit = CDSK::DiffDynSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	s = 10,
	r = 28,
	b = 2.667,
	
	#################
	## Constructor ##
	#################
	
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

