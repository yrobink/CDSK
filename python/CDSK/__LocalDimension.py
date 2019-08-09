# -*- coding: utf-8 -*-


###############
## Libraries ##
###############

import numpy as np
import sys
import scipy.stats as sc
import scipy.spatial.distance as ssd
import sklearn.metrics.pairwise as skmp
import multiprocessing as mp
import SDFC as sd


###########
## Class ##
###########

def _genpareto_fit( distThX , pareto_fit ):
	if pareto_fit == "SDFC":
		law = sd.GPDLaw()
		law.fit( distThX , loc = 0 )
		return law.coef_[0]
	elif pareto_fit == "scipy":
		return sc.genpareto.fit( distThX , floc = 0 )[2]
	else:
		return np.mean( distThX )

## Nicholas Moloney code, original name is extremal_sueveges
def _theta_sueveges_fit( iThreshold , q ):
	Nc = np.count_nonzero( (iThreshold[1:] - iThreshold[:-1] - 1) > 0 )
	N = iThreshold.size - 1
	tmp = ( 1.0 - q ) * ( iThreshold[-1] - iThreshold[0] )
	return ( tmp + N + Nc - np.sqrt( np.power( tmp + N + Nc , 2. ) - 8. * Nc * tmp ) ) / ( 2. * tmp )

def _theta_ferro_fit( iThreshold ):
	Ti = iThreshold[1:] - iThreshold[:-1]
	res = None
	if np.max(Ti) > 2:
		res = 2 * ( np.sum(Ti - 1)**2 ) / ( (Ti.size-1) * np.sum( (Ti-1) * (Ti-2) ) )
	else:
		res = 2 * ( np.sum(Ti)**2 ) / ( (Ti.size-1) * np.sum( Ti**2 ) )
	res = min( 1 , res )
	return res

def _theta_fit( iThreshold , q , theta_fit ):
	if theta_fit == "sueveges":
		return _theta_sueveges_fit( iThreshold , q )
	else:
		return _theta_ferro_fit( iThreshold )

def localDimension_fit( queue , distXY , q , pareto_fit , theta_fit ):
	threshold = np.percentile( distXY , 100 * q , axis = 1 )
	size = threshold.size

	localDim = np.zeros_like( threshold )
	theta = np.zeros_like( threshold )
	
	for i in range(size):
		iThreshold = np.array( np.argwhere( distXY[i,:] > threshold[i] ) ).ravel()
		localDim[i] = 1. / _genpareto_fit( distXY[i,iThreshold] - threshold[i] , pareto_fit )
		theta[i] = _theta_fit( iThreshold , q , theta_fit )
	
	queue[0].put( localDim )
	queue[1].put( theta )


def localDimension( X , Y = None , metric = "euclidean" , q = 0.98 , n_jobs = 1 , pareto_fit = "SDFC" , theta_fit = "ferro" , distXY = None ):
	"""
		CDSK.localDimension
		===================

		Description
		-----------
		Fit a dataset to find its local dimension and persistence index
		
		Parameters
		----------
		X          : np.array[ shape = (n_sample,dimension) ]
			Dataset to fit
		Y          : np.array[ shape = (n_sample_2,dimension) ] or None
			Reference to estimate local dimension. If Y is None, Y = X
		metric     : str or callable = "euclidean"
			Metric used between sample of X and Y, see sklearn.metrics.pairwise.pairwise_distances
		q          : float = 0.98
			Quantile used to find the threshold for generalized pareto distribution
		n_jobs     : int = 1
			Number of CPU available.
		pareto_fit : str = "mean"
			Method to fit the scale of generalized pareto law. Options are : "SDFC" (scale estimation), "scipy" (scale estimation, slower) and "mean" (assume shape = 0).
		theta_fit  : str = "ferro"
			Method to fit the theta. "ferro" or "sueveges".
		distXY     : None or np.array[ shape = (n_sample,n_sample) ]
			Pairwise distance between X and Y, if None, sklearn.metrics.pairwise.pairwise_distances is called
		
		Returns
		-------
		LocalDim   : np.array[ shape = (n_sample) ]
			Local dimension of elements of X
		theta      : np.array[ shape = (n_sample) ]
			Persistence (also called extremal index) of elements of X
	"""

	## Distances
	if distXY is None:
		distXY = skmp.pairwise_distances( X , Y = Y , metric = metric , n_jobs = n_jobs )
	distXY[distXY == 0] = sys.float_info.max
	distXY = - np.log( distXY )
	
	## Parameters for parallel fit
	queue = []
	size,_ = distXY.shape
	idx = np.array_split( range(size) , n_jobs )

	## Fit
	if n_jobs > 1:
		lTh = []
		for i in range(n_jobs):
			queue.append( [mp.Queue(),mp.Queue()] )
			lTh.append( mp.Process( target = localDimension_fit , args = ( queue[-1] , distXY[idx[i]] , q , pareto_fit , theta_fit ) ) )
			lTh[-1].start()
		
		for th in lTh:
			th.join()
	else:
		queue.append( [mp.Queue(),mp.Queue()] )
		localDimension_fit( queue[0] , distXY , q , pareto_fit , theta_fit )

	## Extract results from threads
	localDim = np.zeros(size)
	theta = np.zeros(size)
	for i in range(n_jobs):
		localDim[idx[i]] = queue[i][0].get()
		theta[idx[i]] = queue[i][1].get()
	
	return localDim,theta

