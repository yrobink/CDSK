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
from .__DynamicalSystem import DynamicalSystem


###########
## Class ##
###########

class DiscDynSyst(DynamicalSystem):
	"""
	CDSK.DiscDynSyst
	================
	
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
		

