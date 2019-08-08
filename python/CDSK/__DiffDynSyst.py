# -*- coding: utf-8 -*-



###############
## Libraries ##
###############

import numpy           as np
import scipy.integrate as sci

from .__DynamicalSystem import DynamicalSystem


###########
## Class ##
###########

class DiffDynSyst(DynamicalSystem):
	"""
		CDSK.DiffDynSyst
		================
		
		Description
		-----------
		Abstract base class to construct a continuous dynamical system.
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
		self._edo_solver = "scipy"
	
	
	def _solver( self , X0 , time ):
		if self._edo_solver == "RK4":
			solution = np.zeros( (time.size,X0.size) )
			solution[0,:] = X0
			sTime = time.size - 1
			for i in range(sTime):
				h = time[i+1] - time[i]
				k1 = self._equation( solution[i,:] , time[i] )
				k2 = self._equation( solution[i,:] + h * k1 / 2. , time[i] + h / 2. )
				k3 = self._equation( solution[i,:] + h * k2 / 2. , time[i] + h / 2. )
				k4 = self._equation( solution[i,:] + h * k3 , time[i+1] )
				solution[i+1,:] = solution[i,:] + h * ( k1 + 2 * k2 + 2 * k3 + k4 ) / 6.
			return solution
		elif self._edo_solver == "Euler":
			solution = np.zeros( (time.size,X0.size) )
			solution[0,:] = X0
			sTime = time.size - 1
			for i in range(sTime):
				h = time[i+1] - time[i]
				solution[i+1,:] = solution[i,:] + h * self._equation( solution[i,:] , time[i] )
			return solution
		else:
			return sci.odeint( self , X0 , time )
		

