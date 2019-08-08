# -*- coding: utf-8 -*-


###############
## Libraries ##
###############

import numpy as np
from .__DiffDynSyst import DiffDynSyst


###########
## Class ##
###########

class Lorenz84(DiffDynSyst):
	"""
		CDSK.Lorenz84
		=============

		Description
		-----------
		Lorenz 84 Model, as described in [1] (model), [2] (seasonal cycle) and [3] (climate change forcing).

		[1] E. N. Lorenz, "Irregularity : a Fundamental property of the atmosphere", In Tellus A, vol. 36, n. 2, p. 98-110 (1984)
		[2] E. N. Lorenz, "Can chaos and intransitivity lead to interannual variability?", In Tellus A, vol.42, n. 3, p. 378-389 (1990)
		[3] G. Drotos and al, “Probabilistic concepts in a changing climate : a snapshot attractor picture”. In : Jour. Clim., vol. 28, n. 8, p. 3275–3288 (2015)

	"""
	
	class TimeForcing:##{{{
		"""
			CDSK.Lorenz84.TimeForcing
			=========================
			
			Description
			-----------
			Time forcing of Lorenz84 model
		"""
		def constant( t ):
			"""
				Description
				-----------
				Constant forcing, fixed at value 6.
				
				Parameters
				----------
				t    : float
					Time
			"""
			return 6
		
		def cyclic( t ):
			"""
				Description
				-----------
				Seasonnal forcing, described in [1], a "year" is 73 units of time
				
				[1] E. N. Lorenz, "Can chaos and intransitivity lead to interannual variability?", In Tellus A, vol.42, n. 3, p. 378-389 (1990)
				
				Parameters
				----------
				t    : float
					Time
			"""
			return 9.5 + 2. * np.sin( t * 2. * np.pi / 73. )
		
		def linear( t , tcc = 100 * 73 ):
			"""
				Description
				-----------
				Linear forcing, described in [1], starting at time tcc
				
				[1] G. Drotos and al, “Probabilistic concepts in a changing climate : a snapshot attractor picture”. In : Jour. Clim., vol. 28, n. 8, p. 3275–3288 (2015)
				
				Parameters
				----------
				t    : float
					Time
				tcc  : float
					Starting time of forcing, default is 100 * 73 units of time (year 100)
				"""
			return 0 if t < tcc else -2 * ( t - tcc ) / tcc
	##}}}
	
	def __init__( self , a = 0.25 , b = 4. , G = 1. , F = None , size = 1 ):
		"""
			Parameters
			----------
			a    : float
				default = 0.25
			b    : float
				default = 4.
			G    : float
				default = 1.
			F    : callable or string
			   The time forcing function. default = CDSK.Lorenz84.TimeForcing.constant
					If type(F) == str:
					=> "cyclic"        : seasonal cycle ( == Lorenz84.TimeForcing.cyclic )
					=> "linear"        : linear forcing ( == Lorenz84.TimeForcing.linear )
					=> "cyclic-linear" : seasonal cycle + linear forcing ( == Lorenz84.TimeForcing.cyclic + Lorenz84.TimeForcing.linear )
			size : int
			   Numbers of orbits must be computed
			
			Fix initializations
			-------------------
			dim    : Initialized at 3
			bounds : Initialized at np.array( [ [ -1 , -3 , -3 ] , [ 3 , 3 , 3 ] ] )
		"""
		DiffDynSyst.__init__( self , 3 , size , np.array( [ [ -1 , -3 , -3 ] , [ 3 , 3 , 3 ] ] ) )
		self.a = a
		self.b = b
		self.G = G
		self.F = None
		if callable(F):
			self.F = F
		elif type(F) == str:
			if F == "cyclic":
				self.F = self.TimeForcing.cyclic
			elif F == "linear":
				self.F = self.TimeForcing.linear
			elif F == "cyclic-linear":
				self.F = lambda t : self.TimeForcing.cyclic(t) + self.TimeForcing.linear(t)
		else:
			self.F = self.TimeForcing.constant
	
	
	def _equation( self , X , t ):
		dX = np.zeros(X.shape)
		dX[self._i[0]] = - X[self._i[1]]**2 - X[self._i[2]]**2 - self.a * X[self._i[0]] + self.a * self.F(t)
		dX[self._i[1]] = X[self._i[0]] * X[self._i[1]] - self.b * X[self._i[0]] * X[self._i[2]] - X[self._i[1]] + self.G
		dX[self._i[2]] = X[self._i[0]] * X[self._i[2]] + self.b * X[self._i[0]] * X[self._i[1]] - X[self._i[2]]
		return dX


