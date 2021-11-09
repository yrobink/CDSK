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
	Lorenz 84 Model, as described in [1] (model), [2] (seasonal cycle) and [3]
	(climate change forcing).
	
	[1] E. N. Lorenz, "Irregularity : a Fundamental property of the atmosphere",
	    In Tellus A, vol. 36, n. 2, p. 98-110 (1984)
	[2] E. N. Lorenz, "Can chaos and intransitivity lead to interannual
	    variability?", In Tellus A, vol.42, n. 3, p. 378-389 (1990)
	[3] G. Drotos and al, “Probabilistic concepts in a changing climate : a
	    snapshot attractor picture”. In : Jour. Clim., vol. 28, n. 8,
	    p. 3275–3288 (2015)
	
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
			
			[1] G. Drotos and al, “Probabilistic concepts in a changing climate:
			    a snapshot attractor picture”. In : Jour. Clim., vol. 28, n. 8,
			    p. 3275–3288 (2015)
			
			Parameters
			----------
			t    : float
				Time
			tcc  : float
				Starting time of forcing, default is 100 * 73 units of time
				(year 100)
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
		dX[::3]  = - X[1::3]**2 - X[2::3]**2 - self.a * X[::3] + self.a * self.F(t)
		dX[1::3] = X[::3] * X[1::3] - self.b * X[::3] * X[2::3] - X[1::3] + self.G
		dX[2::3] = X[::3] * X[2::3] + self.b * X[::3] * X[1::3] - X[2::3]
		return dX


