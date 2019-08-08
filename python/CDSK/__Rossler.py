# -*- coding: utf-8 -*-


###############
## Libraries ##
###############

import numpy as np
from .__DiffDynSyst import DiffDynSyst

###########
## Class ##
###########

class Rossler(DiffDynSyst):
	"""
		Apyga.dynamic.continuous.Rossler
		================================

		Description
		-----------
		Rossler Model, as described in [1].
		
		[1] O. E. Rossler, « An equation for continuous chaos », Phys. Rev. Lett. A, vol. 57, np5, 1976, p. 397-398
	"""

	def __init__( self , a = 0.1 , b = 0.1 , c = 14 , size = 1 ):
		"""
			Parameters
			----------
			a    : float
			   Default = 0.1
			b    : float
			   Default = 0.1
			c    : float
			   Default = 14
			size : int
			   Numbers of orbits must be computed

			Fix initializations
			-------------------
			dim    : Initialized at 3
			bounds : Initialized at np.array([ [-20,-20,0] [20,20,35] ])
		"""
		DiffDynSyst.__init__( self , 3 , size , np.array( [ [-20,-20,0] , [20,20,35] ] ) )
		self.a = a
		self.b = b
		self.c = c

	def _equation( self , X , t ):
		dX = np.zeros(X.shape)
		dX[self._i[0]] = - X[self._i[1]] - X[self._i[2]]
		dX[self._i[1]] = X[self._i[0]] + self.a * X[self._i[1]]
		dX[self._i[2]] = self.b + X[self._i[2]] * ( X[self._i[0]] - self.c )
		return dX



