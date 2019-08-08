# -*- coding: utf-8 -*-


###############
## Libraries ##
###############

import numpy as np
from .__DynamicalSystem import DynamicalSystem


###########
## Class ##
###########

class DiscDynSyst(DynamicalSystem):
	"""
		Apyga.dynamic.discrete.DiscDynSyst
		==================================

		Description
		-----------
		Abstract base class to construct a discrete dynamical system.
		This class CAN NOT BE USED directly. It requires to be derived.
	"""

	def __init__( self , dim = 0 , size = 0 , bounds = None ):
		"""
			Initialisation of the continuous dynamical system

			Parameters
			----------
			dim    : int
			   Dimension of the phase space.
			size   : int
			   Numbers of initial condition simultaneously computed by the dynamical system
			bounds : np.array[ shape = (2,dim) ]
			   Bounds of a box in phase space where initial condition can be drawn.
					=> bounds[0,:] is the lower bound
					=> bounds[1,:] is the upper bound
			
			Attributes
			----------
			dim    : int
			   Dimension of the phase space.
			size   : int
			   Numbers of initial condition simultaneously computed by the dynamical system
			bounds : np.array[ shape = (2,dim) ]
			   Bounds of a box in phase space where initial condition can be drawn.
		"""
		DynamicalSystem.__init__( self , dim , size , bounds )
		
	def _solver( self , X0 , time ):
		n = int(time)
		result = np.zeros( (n,X0.size) )
		result[0,:] = X0
		for i in range(1,n,1):
			result[i,:] = self._equation( result[i-1,:] , None )
		return result
		

