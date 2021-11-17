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

class Mandelbrot:
	"""
	CDSK.fractal.Mandelbrot
	=======================
	
	Class to generate the Mandelbrot fractal set
	
	"""
	def __init__( self , x = None , y = None , classic_set = None , nx = None , ny = None , maxit = 200 ):
		"""
		Parameters
		----------
		
		Warning : (x,y) or (classic_set,nx [,ny]) must be set!
		
		x : np.array or None
			x points to estimate the Mandelbrot set
		y : np.array or None
			y points to estimate the Mandelbrot set
		classic_set : str or None
			Pre defined set, available values are "set0", "set1", "set2", "set3" and "set4"
		nx          : int or None
			If classic_set is not None, numbers of points in x-axis
		ny          : int or None
			If classic_set is not None, numbers of points in y-axis
		maxit       : int
			Max number of iterations to estimate the Mandelbrot set. Default is 200
		"""
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
		c = np.array( [ complex( self._x[i] , self._y[j] ) for i,j in itt.product(range(self._nx),range(self._ny)) ] )
		z = np.zeros( (self._nx * self._ny) )
		idx = np.zeros_like(z)
		
		for _ in range(self._maxit):
			z = z**2 + c
			idx[ np.abs(z) < 2 ] += 1
		
		self.ratio = idx.reshape( (self._nx,self._ny) ) / self._maxit
	
	def _set_set0( self ):
		self._x = np.linspace( -0.62 , -0.42 , self._nx )
		self._y = np.linspace( -0.7  , -0.5  , self._ny )
	
	def _set_set1( self ):
		self._x = np.linspace( -0.60  , -0.595 , self._nx )
		self._y = np.linspace( -0.665 , -0.66  , self._ny )
	
	def _set_set2( self ):
		self._x = np.linspace( -0.57  , -0.56  , self._nx )
		self._y = np.linspace( -0.646 , -0.636 , self._ny )
	
	def _set_set3( self ):
		self._x = np.linspace( -0.5663 , -0.5656 , self._nx )
		self._y = np.linspace( -0.6394 , -0.6387 , self._ny )
	
	def _set_set4( self ):
		self._x = np.linspace( -1.5 , -0.5 , self._nx )
		self._y = np.linspace( -0.5 , 0.5  , self._ny )

