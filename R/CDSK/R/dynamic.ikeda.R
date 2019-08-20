

###############
## Libraries ##
###############


###############
## Functions ##
###############

## Ikeda {{{

#' Ikeda
#'
#' Ikeda dynamical system
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param R [float]
#'		  Default = 1
#' @param C1 [float]
#'        Default = 0.4
#' @param C2 [float]
#'        Default = 0.9
#' @param C3 [float]
#'        Default = 6.
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
#'   \item{\code{new(s,r,b,size)}}{This method is used to create object of this class with \code{Ikeda}}
#'   \item{\code{orbit(t,X0)}}{Compute the orbit along t starting at X0. If X0 is NULL, it is randomly drawn by randomIC()}
#'   \item{\code{randomIC()}}{Return a random initial condition}
#' }
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
	
	R = 1.,
	C1 = 0.4,
	C2 = 0.9,
	C3 = 6.,
	
	#################
	## Constructor ##
	#################
	
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

