
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
		return( as.vector( base::apply( self$bounds , 1 , function(X) { return( stats::runif( n = self$size , min = X[1] , max = X[2] ) ) } ) ) )
	},
	
	orbit = function( t , X0 = NULL )
	{
		X0 = if( is.null(X0) ) self$randomIC() else X0
		X = private$solver( t , X0 )
		print(dim(X))
		if( self$size > 1 )
		{
			l = dim(X)[1]
			X = base::array( X , base::c( dim(X)[1] , self$dim , self$size ) )
#			X = base::aperm( X , 
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

## Local Dimension {{{

#' Local Dimension and persistence
#'
#' Compute the local dimension and the persistence of a dataset
#'
#' @param X [matrix] A first matrix (time in column, variables in row). Point where you want the local dimension and persistance
#'
#' @param Y [matrix] A second matrix (time in column, variables in row). Point to estimate ld and theta. If Y = NULL, X is used
#'
#' @param q [float] Threshold, default = 0.98
#'
#' @param distXY [matrix] -log of pairwise distances between X and Y. If NULL, computed with CDSK::pairwise_distances( X , Y , metric = "logeuclidean" )
#'
#' @param gpd_fit [NULL or function] Function which fit the scale parameter of a gpd distribution, take a vector containing in fist index the threshold, and other values the dataset, and return the scale. If NULL, mean inverse is used.
#'
#' @return ld,theta [list] list containing local dimension and theta
#'
#' @examples
#' l63 = CDSK::Lorenz63$new()
#' t = base::seq( 0 , 100 , 0.005 )
#' Y = l63$orbit(t)
#' X = Y[,sample(1:length(t),1000)]
#' ldt = CDSK::localDimension( X , Y )
#' print(base::mean(ldt$ld))
#'
#' @export
localDimension = function( X , Y = NULL , q = 0.98 , distXY = NULL , gpd_fit = NULL )
{
	## Fit function for theta
	theta_ferro = function( Z )
	{
		iThreshold = which( Z[-1] > Z[1] )
		l = length(iThreshold)
		Ti = base::diff(iThreshold)
		res = 2
		if( base::max(Ti) > 2 )
		{
			res = 2 * ( base::sum(Ti - 1)^2 ) / ( (l-1) * base::sum( (Ti-1) * (Ti-2) ) )
		}
		else
		{
			res = 2 * base::sum(Ti)^2 / ( (l-1) * base::sum( Ti^2 ) )
		}
		res = base::min( 1 , res )
		return( res )
	}
	
	## Fit function for local dim
	gpdfit_mean = function(Z)
	{
		return( 1. / base::mean( Z[-1][Z[-1] > Z[1]] - Z[1] ) )
	}
	
	if( is.null(gpd_fit) )
	{
		gpd_fit = gpdfit_mean
	}
	
	## Pairwise distances
	if( is.null(distXY) )
	{
		distXY = - CDSK::pairwise_distances( X , Y , metric = "logeuclidean" )
		distXY[which( is.infinite(distXY) )] = - Inf
	}
	
	## Threshold
	thresholds = base::apply( distXY , 1 , quantile , probs = q )
	
	## Fit localdim
	ld = base::apply( base::cbind( thresholds , distXY ) , 1 , gpd_fit )
#	ld = base::apply( base::cbind( thresholds , distXY ) , 1 , function( Z ) { return(1. / base::mean( Z[-1][Z[-1] > Z[1]] - Z[1] ) ) } )
	
	## Fit theta
	theta = base::apply( base::cbind( thresholds , distXY ) , 1 , theta_ferro )
	
	return( list( ld = ld , theta = theta ) )
}
##}}}


