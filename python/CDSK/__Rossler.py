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

import numpy as np
from .__DiffDynSyst import DiffDynSyst

###########
## Class ##
###########

class Rossler(DiffDynSyst):
	"""
	CDSK.Rossler
	============
	
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
		dX[::3] = - X[1::3] - X[2::3]
		dX[1::3] = X[::3] + self.a * X[1::3]
		dX[2::3] = self.b + X[2::3] * ( X[::3] - self.c )
		return dX



