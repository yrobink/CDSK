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
	
