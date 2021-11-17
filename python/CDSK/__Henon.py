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
		Xnext[::2]  = X[1::2] + 1 - self.a * np.power( X[::2] , 2 )
		Xnext[1::2] = self.b * X[::2]
		return Xnext
	
