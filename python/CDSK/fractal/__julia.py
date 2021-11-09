# -*- coding: utf-8 -*-

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the CDSK (Chaotic        ##
## Dynamical System Kit) library. This library makes it possible                ##
## to generate classic (continuous and discrete) attractors, generate the       ##
## Mandelbrot and Julia set, and fit the local dimension.                       ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## CDSK (Chaotic Dynamical System Kit). Cette librairie permet de générer les   ##
## attracteurs classiques (discret comme continue), de générer l'ensemble de    ##
## Julia et de Mandelbrot et d'estimer les dimensions locales.                  ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################



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
	
	def plot( self , show = False ):
		"""
		Plot the Julia set (after run).
		
		Parameters
		----------
		
		show : bool
			Call plt.show() or no. Default is False
		"""
		X,Y = np.meshgrid( self._x , self._y )
		sx = min( max( self._ny / 100 , 5 ) , 10 )
		sy = min( max( self._nx / 100 , 5 ) , 10 )
		
		fig = plt.figure( figsize = (sx,sy) )
		ax  = fig.add_subplot( 1 , 1 , 1 )
		ax.pcolormesh( X , Y , self.ratio.T , cmap = plt.cm.hot_r )
		fig.set_tight_layout(True)
		if show: plt.show()
		return fig
	
	def export( self , ofile ):
		"""
		Export an image of Julia set.
		
		Parameters
		----------
		ofile : str
			Name of file
		"""
		mpl.image.imsave( ofile , np.rot90(self.ratio) , cmap = plt.cm.hot_r )
	
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


