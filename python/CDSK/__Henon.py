# -*- coding: utf-8 -*-


###############
## Libraries ##
###############

import numpy as np
from .__DiscDynSyst import DiscDynSyst

###########
## Class ##
###########

class Henon(DiscDynSyst):
	"""
		CDSK.Henon
		==========

		Description
		-----------
		Henon Attractor
		
	"""
	def __init__( self , a = 1.4 , b = 0.3 , size = 1 ):
		"""
			Description
			-----------
			Initialisation of parameters of the model and the CDSK.DynamicalSystem class

			Parameters
			----------
			a    : float
			   Default = 1.4
			b    : float
			   Default = 0.3
			size : int
			   Numbers of orbits must be computed

			Fix initializations
			-------------------
			dim    : Initialized at 2
			bounds : Initialized at np.array([ [0,0] [0.5,0.5] ])
		"""
		DiscDynSyst.__init__( self , 2 , size , np.array([ [0,0] , [0.5,0.5] ]) )
		self.a = a
		self.b = b
	
	def _equation( self , X , t = None ):
		Xnext = np.zeros_like(X)
		Xnext[self._i[0]] = X[self._i[1]] + 1 - self.a * np.power( X[self._i[0]] , 2 )
		Xnext[self._i[1]] = self.b * X[self._i[0]]
		return Xnext
	
