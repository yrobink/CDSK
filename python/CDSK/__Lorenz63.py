# -*- coding: utf-8 -*-


###############
## Libraries ##
###############

import numpy as np
from .__DiffDynSyst import DiffDynSyst


###########
## Class ##
###########

class Lorenz63(DiffDynSyst):
	"""
		CDSK.Lorenz63
		=============

		Description
		-----------
		Lorenz 63 Model, as described in [1].

		[1] E. N. Lorenz, « Deterministic nonperiodic flow », J. Atmos. Sci., vol. 20, no 2, 1963, p. 130-141
	"""
	def __init__( self , s = 10 , r = 28 , b = 2.667 , size = 1 ):
		"""
			Description
			-----------
			Initialisation of parameters of the model and the CDSK.DynamicalSystem class

			Parameters
			----------
			s    : float
			   Number of Prandtl, default = 10
			r    : float
			   Number of Rayleigh
			b    : float, default = 28
			   Ratio of critical values, default = 2.667
			size : int
			   Numbers of orbits must be computed

			Fix initializations
			-------------------
			dim    : Initialized at 3
			bounds : Initialized at np.array([ [-20,-20,0] [20,20,40] ])
		"""
		DiffDynSyst.__init__( self , 3 , size , np.array([ [-20,-20,0] , [20,20,40] ]) )
		self.s = s
		self.r = r
		self.b = b
	
	def _equation( self , X , t ):
		dX = np.zeros_like(X)
		dX[self._i[0]] = self.s * ( X[self._i[1]] - X[self._i[0]] )
		dX[self._i[1]] = self.r * X[self._i[0]] - X[self._i[1]] - X[self._i[0]] * X[self._i[2]] 
		dX[self._i[2]] = X[self._i[0]] * X[self._i[1]] - self.b * X[self._i[2]]
		return dX
	
	def jacobian( self , X , t ):
		jac = np.array( [ [ - self.s      , self.s , 0.       ] ,
				            [ self.r - X[2] , - 1.   , - X[0]   ] ,
								[ X[1]          , X[0]   , - self.b ] ] )
		return jac
