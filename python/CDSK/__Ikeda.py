# -*- coding: utf-8 -*-


###############
## Libraries ##
###############

import numpy as np
from .__DiscDynSyst import DiscDynSyst

###########
## Class ##
###########

class Ikeda(DiscDynSyst):
	"""
		Apyga.dynamic.discrete.Ikeda
		============================

		Description
		-----------
		Ikeda attractor.
	"""
	def __init__( self , R = 1. , C1 = 0.4 , C2 = 0.9 , C3 = 6. , size = 1 ):
		"""
			Description
			-----------
			Initialisation of parameters of the model and the Apyga.Dynamic.DynamicalSystem class

			Parameters
			----------
			R    : float
			   Default = 1.
			C1   : float
			   Default = 0.4
			C2   : float
			   Default = 0.9
			C3   : float
			   Default = 6.
			size : int
			   Numbers of orbits must be computed

			Fix initializations
			-------------------
			dim    : Initialized at 2
			bounds : Initialized at np.array([ [0,-1] [1,0] ])
		"""
		DiscDynSyst.__init__( self , 2 , size , np.array([ [ 0. , -1. ] , [ 1. , 0. ] ]) )
		self.R = R
		self.C1 = C1
		self.C2 = C2
		self.C3 = C3
	
	def _equation( self , X , t = None ):
		Xnext = np.zeros_like(X)
		tau = self.C1 - self.C3 / ( 1. + np.power( X[self._i[0]] , 2 ) + np.power( X[self._i[1]] , 2 ) )
		Xnext[self._i[0]] = self.R + self.C2 * ( X[self._i[0]] * np.cos(tau) - X[self._i[1]] * np.sin(tau) )
		Xnext[self._i[1]] = self.C2 * ( X[self._i[1]] * np.cos(tau) + X[self._i[0]] * np.sin(tau) )
		return Xnext
	
