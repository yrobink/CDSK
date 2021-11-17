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

class Ikeda(DiscDynSyst):
	"""
	CDSK.Ikeda
	==========
	
	Description
	-----------
	Ikeda attractor.
	"""
	def __init__( self , R = 1. , C1 = 0.4 , C2 = 0.9 , C3 = 6. , size = 1 ):
		"""
		Description
		-----------
		Initialisation of parameters of the model and the CDSK.DynamicalSystem class
		
		Parameters
		----------
		R    : float
		   Default = 1.
		C1   : float
		   Default = 0.4
		C2   : float
		   Default = 0.9
		C3   : float
		   Default = 6.
		size : int
		   Numbers of orbits must be computed
		
		Fix initializations
		-------------------
		dim    : Initialized at 2
		bounds : Initialized at np.array([ [0,-1] [1,0] ])
		"""
		DiscDynSyst.__init__( self , 2 , size , np.array([ [ 0. , -1. ] , [ 1. , 0. ] ]) )
		self.R = R
		self.C1 = C1
		self.C2 = C2
		self.C3 = C3
	
	def _equation( self , X , t = None ):
		Xnext = np.zeros_like(X)
		tau = self.C1 - self.C3 / ( 1. + X[::2]**2 + X[1::2]**2 )
		Xnext[::2]  = self.R + self.C2 * ( X[::2] * np.cos(tau) - X[1::2] * np.sin(tau) )
		Xnext[1::2] = self.C2 * ( X[1::2] * np.cos(tau) + X[::2] * np.sin(tau) )
		return Xnext
	
