
###############
## Libraries ##
###############


###############
## Functions ##
###############

## Local Dimension {{{

#' Local Dimension and persistence
#'
#' Compute the local dimension and the persistence of a dataset
#'
#' @param X [matrix] A first matrix (time in row, variables in columns). Point where you want the local dimension and persistance
#'
#' @param Y [matrix] A second matrix (time in row, variables in columns). Point to estimate ld and theta. If Y = NULL, X is used
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


