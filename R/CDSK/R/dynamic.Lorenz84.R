

###############
## Libraries ##
###############


###############
## Functions ##
###############

## Lorenz84TimeForcing {{{

#' Lorenz84TimeForcing
#'
#' Lorenz84 time forcing, can be constant, cyclic (period = 73), linear, or cyclic and linear.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param tcc [float]
#'        Time of Climate Change, default = 100 * 73 (100 years)
#' @param t [vector]
#'        Time to evaluate the forcing
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(tcc)}}{This method is used to create object of this class with \code{Lorenz84TimeForcing}}
#'   \item{\code{constant(t)}}{Constant time forcing, fixed at 6}
#'   \item{\code{cyclic(t)}}{Cyclic time forcing, with period fixed at 73, varying between 7.5 and 11.5}
#'   \item{\code{linear(t)}}{Linear forcing, 0 before tcc, and decreasing linearly after tcc.}
#' }
#' @examples
#' ## No example, used by Lorenz84 model
#'
#' @export
Lorenz84TimeForcing = R6::R6Class( "Lorenz84TimeForcing" ,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	tcc = 0,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( tcc = 100 * 73 )
	{
		self$tcc = tcc
	},
	
	
	#############
	## Methods ##
	#############
	
	constant = function(t)
	{
		return(6)
	},
	
	cyclic = function(t)
	{
		return( 9.5 + 2 * base::sin( 2 * base::pi / 73 * t ) )
	},
	
	linear = function( t )
	{
		if( t < self$tcc )
		{
			return(0)
		}
		else
		{
			return( - 2 * ( t - self$tcc ) / self$tcc )
		}
	}
	
	)
)
##}}}

## Lorenz84 {{{

#' Lorenz84
#'
#' Lorenz84 dynamical system
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param a [float]
#'        Default = 0.25
#' @param b [float]
#'        Default = 4.
#' @param G [float]
#'        Default = 1.
#' @param F [callable or string]
#'        Time forcing, if string:
#'        => "constant" use Lorenz84TimeForcing$constant
#'        => "cyclic" use Lorenz84TimeForcing$cyclic
#'        => "linear" use Lorenz84TimeForcing$linear
#'        => "cyclic-linear" use Lorenz84TimeForcing$cyclic + Lorenz84TimeForcing$linear
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
#'   \item{\code{new(a,b,G,F,size)}}{This method is used to create object of this class with \code{Lorenz84}}
#'   \item{\code{orbit(t,X0)}}{Compute the orbit along t starting at X0. If X0 is NULL, it is randomly drawn by randomIC()}
#'   \item{\code{randomIC()}}{Return a random initial condition}
#' }
#' @examples
#' l84 = CDSK::Lorenz84$new( size = 200 , F = "cyclic" )
#' t = base::seq( 0 , 100 , 0.005 )
#' X = l84$orbit(t) 
#' ## X is an array with dim = (length(t),3,200)
#' ## Each X[i,,] is an orbit
#' ## Each X[,,i] is a snapshot
#'
#' @export
Lorenz84 = R6::R6Class( "Lorenz84" , 
	
	inherit = CDSK::DiffDynSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	a = 0.25,
	b = 4.,
	G = 1.,
	F = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( a = 0.25 , b = 4. , G = 1. , F = NULL , size = 1 )
	{
		super$initialize( 3 , size , base::matrix( base::c( -1 , -3 , -3 , 3 , 3 , 3 ) , nrow = 3 , ncol = 2 ) )
		self$a = a
		self$b = b
		self$G = G
		if( is.function(F) )
		{
			self$F = F
		}
		else if( is.character(F) )
		{
			private$forcing = CDSK::Lorenz84TimeForcing$new()
			if( F == "cyclic" )
			{
				self$F = private$forcing$cyclic
			}
			else if( F == "linear" )
			{
				self$F = private$forcing$linear
			}
			else if( F == "cyclic-linear" )
			{
				self$F = function(t) { return( private$forcing$cyclic(t) + private$forcing$linear(t) ) }
			}
		}
		else
		{
			private$forcing = CDSK::Lorenz84TimeForcing$new()
			self$F = private$forcing$constant
		}
	}
	
	
	#############
	## Methods ##
	#############
	
	),
	
	private = list(
	
	###############
	## arguments ##
	###############
	
	forcing = NULL,
	
	
	#############
	## Methods ##
	#############
	
	equation = function( t , X , par = NULL )
	{
		dX = numeric( length(X) )
		dX[private$i[[1]]] = - X[private$i[[2]]]^2 - X[private$i[[3]]]**2 - self$a * X[private$i[[1]]] + self$a * self$F(t)
		dX[private$i[[2]]] = X[private$i[[1]]] * X[private$i[[2]]] - self$b * X[private$i[[1]]] * X[private$i[[3]]] - X[private$i[[2]]] + self$G
		dX[private$i[[3]]] = X[private$i[[1]]] * X[private$i[[3]]] + self$b * X[private$i[[1]]] * X[private$i[[2]]] - X[private$i[[3]]]
		return(list(dX))
	}
	
	)
)
##}}}


