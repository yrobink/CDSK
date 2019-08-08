# -*- coding: utf-8 -*-


###############
## Libraries ##
###############

import numpy as np
from .__DiscDynSyst import DiscDynSyst


###########
## Class ##
###########

class Mira(DiscDynSyst):
	"""
		CDSK.Mira
		=========

		Description
		-----------
		Ikeda attractor.
	"""
	def __init__( self , a = -0.48 , b = 0.93 , size = 1 ):
		"""
			Description
			-----------
			Initialisation of parameters of the model and the CDSK.DynamicalSystem class

			Parameters
			----------
			a    : float
			   Default = -0.48
			b    : float
			   Default = 0.93
			size : int
			   Numbers of orbits must be computed

			Fix initializations
			-------------------
			dim    : Initialized at 2
			bounds : Initialized at np.array([ [ 3.8 , -0.2 ] , [ 4.2 , 0.2 ] ])
		"""
		DiscDynSyst.__init__( self , 2 , size , np.array([ [ 3.8 , -0.2 ] , [ 4.2 , 0.2 ] ]) )
		self.a = a
		self.b = b
	
	def _F( self , x ):
		return self.a * x + 2 * ( 1 - self.a ) * np.power( x , 2 ) / ( 1 + np.power( x , 2 ) )
	
	def _equation( self , X , t = None ):
		Xnext = np.zeros_like(X)
		Xnext[self._i[0]] = self.b * X[self._i[1]] + self._F( X[self._i[0]] )
		Xnext[self._i[1]] = - X[self._i[0]] + self._F( Xnext[self._i[0]] )
		return Xnext
	
