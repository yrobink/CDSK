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
import itertools as itt
import matplotlib as mpl
import matplotlib.pyplot as plt


###############
## Functions ##
###############

class Julia:
	"""
	CDSK.fractal.Julia
	==================
	
	Class to generate the Julia fractal set
	
	"""
	def __init__( self , c = None , x = None , y = None , classic_set = None , nx = None , ny = None , maxit = 200 ):
		"""
		Parameters
		----------
		
		Warning : (x,y) or (classic_set,nx [,ny]) must be set!
		
		x : np.array or None
			x points to estimate the Julia set
		y : np.array or None
			y points to estimate the Julia set
		classic_set : str or None
			Pre defined set, available values are "set0", "set1", "set2", "set3" and "set4"
		nx          : int or None
			If classic_set is not None, numbers of points in x-axis
		ny          : int or None
			If classic_set is not None, numbers of points in y-axis
		maxit       : int
			Max number of iterations to estimate the Julia set. Default is 200
		"""
		self._c           = c
		self._x           = x
		self._y           = y
		self._classic_set = classic_set
		self._nx          = nx
		self._ny          = ny
		self._maxit       = maxit
		self.ratio        = None
		
		try:
			if self._classic_set is not None and (nx is not None or ny is not None):
				
				self._nx = self._nx if self._nx is not None else self._ny
				self._ny = self._ny if self._ny is not None else self._nx
				
				if   self._classic_set == "set0": self._set_set0()
				elif self._classic_set == "set1": self._set_set1()
				elif self._classic_set == "set2": self._set_set2()
				elif self._classic_set == "set3": self._set_set2()
				elif self._classic_set == "set4": self._set_set4()
				
			else:
				self._nx = self._x.size
				self._ny = self._y.size
		except:
			raise NameError("Need (x,y) set or (classic_set,nx [,ny]) set!")
	
	@property
	def x(self):
		return self._x
	
	@property
	def y(self):
		return self._y
	
	def run(self):
		"""
		Run estimation
		"""
		z = np.array( [ complex( self._x[i] , self._y[j] ) for i,j in itt.product(range(self._nx),range(self._ny)) ] )
		idx = np.zeros( (self._nx * self._ny) , dtype = np.float )
		
		for _ in range(self._maxit):
			z = z**2 + self._c
			idx[ np.abs(z) < 2 ] += 1
		
		self.ratio = idx.reshape( (self._nx,self._ny) ) / self._maxit
	
	def _set_set0( self ):
		self._c = complex(0.3,0.5)
		self._x = np.linspace( 0.25 , 0.53 , self._nx )
		self._y = np.linspace( 0.4  , 0.72 , self._ny )
	
	def _set_set1( self ):
		self._c = complex(0.285,0.01)
		self._x = np.linspace( -1   , 1   , self._nx )
		self._y = np.linspace( -1.1 , 1.1 , self._ny )
	
	def _set_set2( self ):
		self._c = complex(0.285,0.013)
		self._x = np.linspace( -1   , 1   , self._nx )
		self._y = np.linspace( -1.1 , 1.1 , self._ny )
	
	def _set_set3( self ):
		self._c = complex( -0.038088 , 0.9754633 )
		self._x = np.linspace( -2e-3 , 2e-3 , self._nx )
		self._y = np.linspace( -2e-3 , 2e-3 , self._ny )
	
	def _set_set4( self ):
		self._c = complex( -0.8 , 0.156 )
		self._x = np.linspace( -1.55 , 1.55 , self._nx )
		self._y = np.linspace( -1    , 1    , self._ny )


