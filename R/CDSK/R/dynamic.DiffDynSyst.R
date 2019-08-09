
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
#' ## X is an array with dim = (length(t),3,200)
#' ## Each X[i,,] is an orbit
#' ## Each X[,,i] is a snapshot
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


