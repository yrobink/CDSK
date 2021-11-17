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


###########
## Class ##
###########

class DynamicalSystem:
	"""
	CDSK.DynamicalSystem
	====================
	
	Description
	-----------
	Abstract base class to construct a dynamical system.
	This class CAN NOT BE USED directly. It requires to be derived.
	"""

	def __init__( self , dim = 0 , size = 0 , bounds = None ):
		"""
		Initialisation of the dynamical system
		
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
		self.dim = dim
		self.size = size
		self.bounds = bounds
	
	def randomIC( self ):
		"""
			Draw uniformly an Initial Condition in a box defined by self.bounds
		"""
		return np.random.uniform( low = self.bounds[0,:] , high = self.bounds[1,:] , size = ( self.size , self.dim ) )
	
	def _equation( self , X , t ):
		pass
	
	def __call__( self , X , t ):
		"""
		Return the value of the equation of the dynamical system at point X (and time t if continuous)
		
		Parameters
		----------
		X   : np.array[ shape = (self.size * self.dim) ]
		   Vector to computes the derivative.
		t   : float
		   Time to computes the derivative. (not used if system is discrete)
		
		Returns
		-------
		out : array_like
		   Derivative of the dynamical system at point X and time t
		"""
		return self._equation( X , t )
	
	def _solver( self , time , X0 ):
		pass
	
	def jacobian( self , X , t ):
		pass

	def orbit( self , time , X0 = None ):
		"""
		Build an orbit of the dynamical system. If self.size > 1, a snapshot attractor [1] is built.
		
		Parameters
		----------
		time : array_like
		   An array of time to estimate the value of dynamical system
		X0   : np.array[ shape = ( self.size , self.dim ) ]
		   Initial condition corresponding to time[0].
				=> If X0 is None, self.randomIC() is called
		
		Returns
		-------
		out  : array_like
		   An array of the orbit. The shape is (len(time),self.size,self.dim)
		
		Ref
		---
		[1] Romeiras, F. J., Grebogi, C., and Ott, E.: Multifractal properties
		    of snapshot attractors of random maps, Phys. Rev A, 41, 784–799,
		    doi:10.1103/PhysRevA.41.784
		"""
		if X0 is None:
			X0 = self.randomIC()
		
		if len(X0.shape) > 1:
			X0 = X0.reshape( (X0.size) )
		
		_orbit = self._solver( X0 , time )
		
		lTime = time if np.isscalar(time) else len(time)
		
		if self.size == 1:
			_orbit = _orbit.reshape( (lTime,self.dim) )
		else:
			_orbit = _orbit.reshape( (lTime,self.size,self.dim) )
		
		return _orbit
	
