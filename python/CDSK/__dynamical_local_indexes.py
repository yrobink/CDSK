# -*- coding: utf-8 -*-

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

###############
## Libraries ##
###############

import sys,os
import itertools as itt
import numpy as np
import scipy.stats as sc
import sklearn.metrics.pairwise as skmp
import multiprocessing as mp
import SDFC as sd


###############
## Functions ##
###############

def _theta_sueveges( idx , ql ):##{{{
	N   = idx.size - 1
	Nc  = np.sum( idx[1:] - idx[:-1] - 1 > 0 )
	T   = ( 1. - ql ) * ( idx[-1] - idx[0] )
	return ( T + N + Nc - np.sqrt( np.power( T + N + Nc , 2. ) - 8. * Nc * T ) ) / ( 2. * T )
##}}}

def _theta_ferro( idx ): ##{{{
	Ti = idx[1:] - idx[:-1]
	if np.max(Ti) > 2:
		res = 2 * ( np.sum(Ti - 1)**2 ) / ( (Ti.size-1) * np.sum( (Ti-1) * (Ti-2) ) )
	else:
		res = 2 * ( np.sum(Ti)**2 ) / ( (Ti.size-1) * np.sum( Ti**2 ) )
	return min( 1 , res )
##}}}

def _gpd_fit( data , ld_fit ):##{{{
	if ld_fit == "scipy":
		shape,_,scale = sc.genpareto.fit( data , floc = 0 )
		return np.array( [scale,shape] )
	else:
		gpd = sd.GPD()
		gpd.fit( data , f_loc = 0 )
		return gpd.coef_
##}}}

def _local_dimension_parallel( queue , dist , q , ld_fit ): ##{{{
	size = dist.shape[0]
	ld   = np.zeros(size)
	shp  = np.zeros(size)
	
	for i in range(size):
		data = dist[i,dist[i,:] > q[i]] - q[i]
		ld[i],shp[i] = _gpd_fit( data , ld_fit )
	
	queue[0].put(ld)
	queue[1].put(shp)
##}}}

def _local_dimension( dist , q , where , ld_fit , n_jobs ):##{{{
	
	## Mean estimator
	if ld_fit not in ["SDFC","scipy"]:
		return 1. / np.nanmean( np.where( where , dist - q , np.nan ) , 1 ),np.zeros(q.size)
	
	## Others
	n_sample  = q.size
	l_idx     = np.array_split( range(n_sample) , n_jobs )
	l_queue   = []
	l_threads = []
	for idx in l_idx:
		l_queue.append( [mp.Queue(),mp.Queue()] )
		l_threads.append( mp.Process( target = _local_dimension_parallel , args = (l_queue[-1],dist[idx,:],q[idx],ld_fit) ) )
		l_threads[-1].start()
	
	for th in l_threads:
		th.join()
	

	ld  = np.zeros(n_sample)
	shp = np.zeros(n_sample)
	
	for idx,queue in zip(l_idx,l_queue):
		ld[idx]  = 1. / queue[0].get()
		shp[idx] = queue[1].get()
	
	return ld,shp
##}}}

def dynamical_local_indexes( X , Y = None , ql = 0.98 , ld_fit = "SDFC" , theta_fit = "sueveges" , n_jobs = os.cpu_count() , **kwargs ):##{{{
	"""
	CDSK.localDimension
	===================
	
	Description
	-----------
	Fit a dataset to find (cross) local dimension, persistence, alpha and shape.
	
	Parameters
	----------
	X          : np.array[ shape = (n_sample,n_features,n_var) ]
		Dataset to fit
	Y          : np.array[ shape = (n_sample_2,n_features,n_var) ] or None
		Reference to estimate. If Y is None, Y = X
	ql         : float = 0.98
		Quantile used to find the threshold for generalized pareto distribution
	ld_fit     : str = "SDFC"
		Method to fit the scale of generalized pareto law. Options are : "SDFC"
		(scale estimation), "scipy" (scale estimation, slower) and "mean"
		(assume shape = 0).
	theta_fit  : str = "sueveges"
		Method to fit the persistence. "ferro" or "sueveges".
	n_jobs     : int = os.cpu_count()
		Number of CPU available.
	metric     : str or callable = "euclidean"
		Metric used between sample of X and Y, see sklearn.metrics.pairwise.pairwise_distances
	cross_metric: callable = lambda x,y : np.sqrt(x**2+y**2)
		Metric used between two variables.
	
	
	Returns
	-------
	ld    : np.array[ shape = (n_sample,n_var,n_var) ]
		Local dimension of elements of X. The vector ld[:,i,j] is the cross-local
		dimension between X[:,:,i] and X[:,:,j]. Only the upper part of the 
		matrix is filled.
	theta : np.array[ shape = (n_sample,n_var,n_var) ]
		Persistence (also called extremal index) of elements of X. The vector
		theta[:,i,j] is the cross persistence between X[:,:,i] and X[:,:,j]. 
		Only the upper part of the matrix is filled.
	alpha : np.array[ shape = (n_sample,n_var,n_var) ]
		Co-recurrences of elements of X. The vector alpha[:,i,j] is the
		co-recurrences between X[:,:,i] and X[:,:,j].  Only the upper part of
		the matrix is filled. See [1].
	shape : np.array[ shape = (n_sample,n_var,n_var) ]
		Shape estimated for the GPD during the fit of local dimension. The
		vector shape[:,i,j] is the shape between X[:,:,i] and X[:,:,j].  Only
		the upper part of the matrix is filled. If ld_fit == "mean", shape is
		always 0.
	
	References
	----------
	[1] Messori. G. and Faranda, D. (2021) Technical note: Characterising and
	comparing different palaeoclimates with dynamical systems theory.
	Clim. Past, 17, 545–563. doi:https://doi.org/10.5194/cp-17-545-2021.
	"""
	
	## Read kwargs
	metric = kwargs.get("metric")
	if metric is None: metric = "euclidean"
	cross_metric = kwargs.get("cross_metric")
	if cross_metric is None: cross_metric = lambda x,y : np.sqrt(x**2+y**2)
	
	## Data in good format
	if X.ndim == 2:
		X = X.reshape( X.shape + (1,) )
	if Y is None:
		Y = X
	elif Y.ndim == 2:
		Y = Y.reshape( Y.shape + (1,) )
	
	n_sampleX,n_features,n_traj = X.shape
	n_sampleY,n_features,n_traj = Y.shape
	
	## Distances
	dist = np.zeros( (n_sampleX,n_sampleY,n_traj,n_traj) )
	for i in range(n_traj):
		dist[:,:,i,i] = skmp.pairwise_distances( X[:,:,i] , Y[:,:,i] , metric = metric , n_jobs = n_jobs )
		dist[:,:,i,i] /= np.linalg.norm( dist[:,:,i,i] , 1 )
	
	for i,j in itt.combinations(range(n_traj),2):
		dist[:,:,i,j] = cross_metric(dist[:,:,i,i],dist[:,:,j,j])
	
	dist[~(dist > 0)] = sys.float_info.max
	dist = -np.log(dist)
	
	## Local dimension
	ld    = np.zeros( (n_sampleX,n_traj,n_traj) ) + np.nan
	shp   = np.zeros( (n_sampleX,n_traj,n_traj) ) + np.nan
	q     = np.zeros( (n_sampleX,n_traj,n_traj) ) + np.nan
	where = np.zeros( (n_sampleX,n_sampleY,n_traj,n_traj) , dtype = bool )
	for i,j in itt.combinations_with_replacement(range(n_traj),2):
		q[:,i,j]             = np.quantile( dist[:,:,i,j] , ql , 1 )
		where[:,:,i,j]       = dist[:,:,i,j] > q[:,i,j].reshape(-1,1)
		ld[:,i,j],shp[:,i,j] = _local_dimension( dist[:,:,i,j] , q[:,i,j].reshape(-1,1) , where[:,:,i,j] , ld_fit , n_jobs )
	
	## Theta
	theta = np.zeros( (n_sampleX,n_traj,n_traj) ) + np.nan
	for i,j in itt.combinations_with_replacement(range(n_traj),2):
		for k in range(n_sampleX):
			idx = np.argwhere(where[k,:,i,j]).squeeze()
			if theta_fit == "ferro":
				theta[k,i,j] = _theta_ferro( idx , ql )
			else:
				theta[k,i,j] = _theta_sueveges( idx , ql )
	
	## Alpha
	alpha = np.zeros( (n_sampleX,n_traj,n_traj) ) + np.nan
	for i,j in itt.combinations_with_replacement(range(n_traj),2):
		alpha[:,i,j] = np.sum( where[:,:,i,i] & where[:,:,j,j] , 1 ) / np.sum(where[:,:,i,i],1)
	
	
	return ld,theta,alpha,shp
##}}}

