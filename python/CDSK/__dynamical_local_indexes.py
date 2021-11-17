# -*- coding: utf-8 -*-

## Copyright(c) 2021 Yoann Robin
## 
## This file is part of CDSK.
## 
## CDSK is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## CDSK is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with CDSK.  If not, see <https://www.gnu.org/licenses/>.

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
	q = 1 - ql
	Ti = idx[1:] - idx[:-1]
	Si = Ti - 1
	Nc = np.sum(Si > 0)
	K  = np.sum( q * Si )
	N  = Ti.size
	return ( K + N + Nc - np.sqrt( (K + N + Nc)**2 - 8 * Nc * K ) ) / (2 * K)

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
	CDSK.dynamical_local_indexes
	============================
	
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
	return_shape: bool = False
		Return the shape estimated by the GPD fit.
	return_dist: bool = False
		Return the pairwise distances.
	return_where: bool = False
		Return the 'where' matrix.
	
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
	
	Optional returns
	----------------
	shape : np.array[ shape = (n_sample,n_var,n_var) ]
		Shape estimated for the GPD during the fit of local dimension. The
		vector shape[:,i,j] is the shape between X[:,:,i] and X[:,:,j].  Only
		the upper part of the matrix is filled. If ld_fit == "mean", shape is
		always 0.
	dist  : np.array[ shape = (n_sample,n_sample_2,n_var,n_var) ]
		Pairwise distances between X and Y.
	where : np.array[ shape = (n_sample,n_sample_2,n_var,n_var) , dtype = bool ]
		True if distance between X/Y is greater than the quantile ql, otherwise
		0.
	
	References
	----------
	[1] Messori. G. and Faranda, D. (2021) Technical note: Characterising and
	comparing different palaeoclimates with dynamical systems theory.
	Clim. Past, 17, 545–563. doi:https://doi.org/10.5194/cp-17-545-2021.
	[2] Süveges, Mária. 2007. Likelihood estimation of the extremal index.
	Extremes, 10.1-2, 41-55, doi:10.1007/s10687-007-0034-2
	"""
	
	## Read kwargs
	metric = kwargs.get("metric")
	if metric is None: metric = "euclidean"
	cross_metric = kwargs.get("cross_metric")
	if cross_metric is None: cross_metric = lambda x,y : np.sqrt(x**2+y**2)
	return_shape = kwargs.get("return_shape")
	return_dist  = kwargs.get("return_dist")
	return_where = kwargs.get("return_where")
	if return_shape is None:
		return_shape = False
	if return_dist is None:
		return_dist = False
	if return_where is None:
		return_where = False
	
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
				theta[k,i,j] = _theta_ferro( idx )
			else:
				theta[k,i,j] = _theta_sueveges( idx , ql )
	
	## Alpha
	alpha = np.zeros( (n_sampleX,n_traj,n_traj) ) + np.nan
	for i,j in itt.combinations_with_replacement(range(n_traj),2):
		alpha[:,i,j] = np.sum( where[:,:,i,i] & where[:,:,j,j] , 1 ) / np.sum(where[:,:,i,i],1)
	
	## Build tuple output
	out = (ld,theta,alpha)
	
	if return_shape:
		out = out + (shp,)
	if return_dist:
		out = out + (dist,)
	if return_where:
		out = out + (where,)
	
	return out
##}}}

