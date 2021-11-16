
##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2021                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the CDSK (Chaotic        ##
## Dynamical System Kit) library. This library makes it possible                ##
## to generate classic (continuous and discrete) attractors, generate the       ##
## Mandelbrot and Julia set, and fit the local dimension.                       ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2021                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## CDSK (Chaotic Dynamical System Kit). Cette librairie permet de générer les   ##
## attracteurs classiques (discret comme continue), de générer l'ensemble de    ##
## Julia et de Mandelbrot et d'estimer les dimensions locales.                  ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################

gpd_fit = function(dNA)##{{{
{
	dNA = dNA[is.finite(dNA)]
	gpd = base::try(ROOPSD::GPD$new()$fit(dNA,loc=0),silent=TRUE)
	if( "try-error" %in% class(gpd) )
	{
		return( base::c( 1. / base::mean(dNA,na.rm=TRUE) , 0 ) )
	}
	return(base::c(gpd$scale,gpd$shape))
}
##}}}

local_dimension = function( dist , q , ld_fit )##{{{
{
	q  = matrix( q , nrow = length(q) )
	dNA = base::apply( dist , 2 , function(x) { x - q } )
	dNA[dNA < 0] = NA
	if( ld_fit == "GPD" )
	{
		ldshp = base::apply( dNA , 1 , gpd_fit )
		ld  = ldshp[1,]
		shp = ldshp[2,]
	}
	else
	{
		ld  = 1. / base::apply( dNA , 1 , base::mean , na.rm = TRUE )
		shp = 0
	}
	
	return( list( ld = ld , shp = shp ) )
}
##}}}

persistence = function( idx , ql , theta_fit )##{{{
{
	Ti = idx[2:length(idx)] - idx[1:(length(idx)-1)]
	N  = length(Ti)
	if( theta_fit == "sueveges" )
	{
		q = 1 - ql
		
		Si = Ti - 1
		Nc = base::sum( Si > 0 )
		K  = base::sum( q * Si )
		th = ( K + N + Nc - base::sqrt( (K + N + Nc)^2 - 8 * Nc * K ) ) / (2 * K)
	}
	else
	{
		if( base::max(Ti) > 2 )
		{
			th = 2 * ( base::sum(Ti - 1)^2 ) / ( (N-1) * base::sum( (Ti-1) * (Ti-2) ) )
		}
		else
		{
			th = 2 * ( base::sum(Ti)^2 ) / ( (N-1) * base::sum( Ti^2 ) )
		}
	}
	
	th[th < 0] = 0
	th[th > 1] = 1
	invisible(th)
}
##}}}

## dynamical_local_indexes ##{{{

#' dynamical_local_indexes
#'
#' Compute the local dimension, persistence and co-reccurrence of a dataset
#'
#' @importFrom pmetric pairwise_distances
#'
#' @param X [array] A first matrix  (n_samples x n_features x n_variables).
#'        Point where you want the local dimension and persistance
#' @param Y [array] A second matrix (n_samples x n_features x n_variables).
#'        Point to estimate ld and theta. If Y = NULL, X is used
#' @param ql [float] Threshold, default = 0.98
#' @param ld_fit [string] Method for compute local_dimension. Currently, only
#'        the mean method is supported, assuming the shape parameters of the GPD
#'        is 0.
#' @param theta_fit [string] Method to fit the persistence. "sueveges" or
#'        "ferro".
#' @param ... Others parameters:
#'        metric [string or function] Metric used for pairwise distances between
#'           X and Y. See the function pmetric::pairwise_distances. Default is
#'           "euclidean".
#'        cross_metric [function] Metric for co-reccurences. Default is the
#'           square of euclidean distance.
#'        return_shape [bool] Return the shape fitted of the GPD distribution.
#'        return_dist [bool] Return the pairwise distances.
#'        return_where [bool] Return the where matrix.
#'
#' @return ld,theta,alpha[,shape,dist,where] [list] list containing local
#'         dimension, persistence, co-recurrences, shape. 'where' is an array
#'         such that where[i,j,v1,v2] is TRUE if the -log of the distance
#'         between X[i,,v1] and Y[j,,v2] is greater than the quantile ql.
#'
#' @examples
#' l63 = CDSK::Lorenz63$new()
#' t = base::seq( 0 , 100 , length = 1000 )
#' O = l63$orbit(t)[501:1000,]
#' X = array( NA , dim = base::c( 500 , 2 , 3 ) )
#' X[,,1] = O[,1:2]
#' X[,,2] = O[,2:3]
#' X[,,3] = O[,base::c(1,3)]
#' dli = CDSK::dynamical_local_indexes(X)
#'
#' @export
dynamical_local_indexes = function( X , Y = NULL , ql = 0.98 , ld_fit = "mean" , theta_fit = "sueveges" , ... )
{
	## Read kwargs
	kwargs = list(...)
	metric = kwargs[["metric"]]
	cross_metric = kwargs[["cross_metric"]]
	return_shape = kwargs[["return_shape"]]
	return_dist  = kwargs[["return_dist"]]
	return_where = kwargs[["return_where"]]
	
	if( is.null(metric) )
		metric = "euclidean"
	if( is.null(cross_metric) )
		cross_metric = function(x,y) { return(x^2+y^2) }
	if( is.null(return_shape) )
		return_shape = FALSE
	if( is.null(return_dist) )
		return_dist = FALSE
	if( is.null(return_where) )
		return_where = FALSE
	
	
	## Read args
	if( is.null(Y) )
		Y = X
	
	n_samplesX = dim(X)[1]
	n_samplesY = dim(Y)[1]
	n_features = dim(X)[2]
	n_traj     = dim(X)[3]
	
	
	## Distances
	dist = array( NA , base::c(n_samplesX,n_samplesY,n_traj,n_traj) )
	
	for( i in 1:n_traj )
	{
		dist[,,i,i] = pmetric::pairwise_distances( X[,,i] , Y[,,i] , metric = metric )
		dist[,,i,i] = dist[,,i,i] / base::apply( dist[,,i,i] , 1 , base::sum )
	}
	
	for( i in 1:n_traj )
	{
		for( j in (i+1):n_traj )
		{
			if( !(j>n_traj) )
				dist[,,i,j] = cross_metric( dist[,,i,i] , dist[,,j,j] )
		}
	}
	dist = -base::log(dist)
	dist[!is.finite(dist)] = NA
	
	## Local dimension
	ld    = array( NA , dim = base::c(n_samplesX,n_traj,n_traj) )
	shp   = array( NA , dim = base::c(n_samplesX,n_traj,n_traj) )
	q     = array( NA , dim = base::c(n_samplesX,n_traj,n_traj) )
	where = array( NA , dim = base::c(n_samplesX,n_samplesY,n_traj,n_traj) )
	
	for( i in 1:n_traj )
	{
		for( j in i:n_traj )
		{
			q[,i,j]      = base::apply( dist[,,i,j] , 1 , stats::quantile , probs = ql , na.rm = TRUE )
			where[,,i,j] = base::apply( dist[,,i,j] , 2 , function(x){ x > q[,i,j] } )
			ldshp        = local_dimension( dist[,,i,j] , q[,i,j] , ld_fit )
			ld[,i,j]     = ldshp$ld
			shp[,i,j]    = ldshp$shp
		}
	}
	
	## Theta
	theta = array( NA , dim = base::c(n_samplesX,n_traj,n_traj) )
	
	for( i in 1:n_traj )
	{
		for( j in i:n_traj )
		{
			for( k in 1:n_samplesX )
			{
				idx = which( where[k,,i,j] )
				theta[k,i,j] = persistence( idx , ql , theta_fit )
			}
		}
	}
	
	## Alpha
	alpha = array( NA , dim = base::c(n_samplesX,n_traj,n_traj) )
	for( i in 1:n_traj )
	{
		for( j in i:n_traj )
		{
			alpha[,i,j] = base::apply( where[,,i,i] & where[,,j,j] , 1 , base::sum ) / base::apply( where[,,i,i] , 1 , base::sum )
		}
	}
	
	## Output
	out = list( ld = ld , theta = theta , alpha = alpha )
	if( return_shape )
		out[["shape"]] = shp
	if( return_dist )
		out[["dist"]] = dist
	if( return_where )
		out[["where"]] = where
	
	invisible(out)
}
##}}}


