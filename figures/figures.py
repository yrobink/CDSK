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



#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

###############
## Libraries ##
###############

import sys,os
import pickle as pk
import multiprocessing as mp

import numpy as np
import pandas as pd
import scipy.stats as sc
import scipy.optimize as sco

import matplotlib as mpl
try:
	import matplotlib.pyplot as plt
except:
	mpl.use("Qt5Agg")
	import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d


import CDSK as ck
import CDSK.fractal as ckf


####################
## Paramètres mpl ##
####################

mpl.rcParams['font.size'] = 15


###############
## Fonctions ##
###############

def test_lorenz63( plot = True ):##{{{
	print( "Test Lorenz63..." , end = "\r" )
	try:
		l63 = ck.Lorenz63()
		X0  = l63.orbit( np.linspace( 0 , 100 , 100 ) )[-1,:]
		X   = l63.orbit( np.linspace( 0 , 100 , 10000 ) , X0 = X0 )
		
		if plot:
			fig = plt.figure( figsize = (7,7) )
			ax = fig.add_subplot( 1 , 1 , 1 , projection = "3d" )
			ax.plot( X[:,0] , X[:,1] , X[:,2] , color = "red" , linestyle = "-" )
			fig.set_tight_layout(True)
			plt.show()
		print( "Test Lorenz63 (Done)  " )
	except:
		print( "Test Lorenz63 (Fail)  " )
##}}}

def test_rossler( plot = True ):##{{{
	print( "Test Rossler..." , end = "\r" )
	try:
		ross = ck.Rossler()
		X0  = ross.orbit( np.linspace( 0 , 100 , 100 ) )[-1,:]
		X   = ross.orbit( np.linspace( 0 , 300 , 10000 ) , X0 = X0 )
		
		if plot:
			fig = plt.figure( figsize = (7,7) )
			ax = fig.add_subplot( 1 , 1 , 1 , projection = "3d" )
			ax.plot( X[:,0] , X[:,1] , X[:,2] , color = "red" , linestyle = "-" )
			fig.set_tight_layout(True)
			plt.show()
		print( "Test Rossler (Done)  " )
	except:
		print( "Test Rossler (Fail)  " )
##}}}

def test_lorenz84( plot = True ):##{{{
	print( "Test Lorenz84..." , end = "\r" )
	try:
		l63 = ck.Lorenz84( size = 1000 , F = "cyclic" )
		X0  = l63.orbit( np.linspace( 0 , 2 * 73 , 50 ) )[-1,:,:]
		X   = l63.orbit( np.linspace( 0 , 5 * 73 , 10000 ) , X0 = X0 )
		t_fall = int( 4    / 5 * 10000 )
		t_wint = int( 4.25 / 5 * 10000 )
		t_spri = int( 4.5  / 5 * 10000 )
		t_summ = int( 4.75 / 5 * 10000 )
		lX = []
		lX.append( X[t_fall,:,:] )
		lX.append( X[t_wint,:,:] )
		lX.append( X[t_spri,:,:] )
		lX.append( X[t_summ,:,:] )
		title = ["Fall","Winter","Spring","Summer"]
		
		if plot:
			fig = plt.figure( figsize = (12,12) )
			for i,X in enumerate(lX):
				ax = fig.add_subplot( 2 , 2 , i + 1 , projection = "3d" )
				ax.plot( X[:,0] , X[:,1] , X[:,2] , color = "red" , linestyle = "" , marker = "." )
				ax.set_title( title[i] , fontsize = 20 )
			fig.set_tight_layout(True)
			plt.show()
		print( "Test Lorenz84 (Done)  " )
	except:
		print( "Test Lorenz84 (Fail)  " )
##}}}

def test_henon( plot = True ):##{{{
	print( "Test Henon..." , end = "\r" )
	try:
		henon = ck.Henon()
		X   = henon.orbit( 10000 , X0 = np.array( [0.5,0.5] ) )
		
		if plot:
			fig = plt.figure( figsize = (7,7) )
			ax = fig.add_subplot( 1 , 1 , 1 )
			ax.plot( X[:,0] , X[:,1] , color = "red" , linestyle = "" , marker = "." )
			fig.set_tight_layout(True)
			plt.show()
		print( "Test Henon (Done)  " )
	except:
		print( "Test Henon (Fail)  " )
##}}}

def test_ikeda( plot = True ):##{{{
	print( "Test Ikeda..." , end = "\r" )
	try:
		ike = ck.Ikeda()
		X   = ike.orbit( 10000 , X0 = np.array( [0.5,0.5] ) )
		
		if plot:
			fig = plt.figure( figsize = (7,7) )
			ax = fig.add_subplot( 1 , 1 , 1 )
			ax.plot( X[:,0] , X[:,1] , color = "red" , linestyle = "" , marker = "." )
			fig.set_tight_layout(True)
			plt.show()
		print( "Test Ikeda (Done)  " )
	except:
		print( "Test Ikeda (Fail)  " )
##}}}

def test_mira( plot = True ):##{{{
	print( "Test Mira..." , end = "\r" )
	try:
		mir = ck.Mira()
		X0  = mir.orbit( 10000 , X0 = np.array( [0.5,0.5] ) )
		mir = ck.Mira( size = 10000 )
		X   = mir.orbit( 1000 , X0 = X0 )[-1,:,:]
		
		if plot:
			fig = plt.figure( figsize = (7,7) )
			ax = fig.add_subplot( 1 , 1 , 1 )
			ax.plot( X[:,0] , X[:,1] , color = "red" , linestyle = "" , marker = "." )
			fig.set_tight_layout(True)
			plt.show()
		print( "Test Mira (Done)  " )
	except:
		print( "Test Mira (Fail)  " )
##}}}


def test_local_dimension( plot = True ):##{{{
	print( "Test local dimension..." , end = "\r" )
	try:
		l63 = ck.Lorenz63()
		X0  = l63.orbit( np.linspace( 0 , 100 , 100 ) )[-1,:]
		X   = l63.orbit( np.linspace( 0 , 100 , 10000 ) , X0 = X0 )
		
		ld,th = ck.localDimension( X , n_jobs = 6 , pareto_fit = "SDFC" )
		vmin_ld,vmax_ld = np.quantile( ld , [0.1,0.9] )
		vmin_th,vmax_th = np.quantile( th , [0.1,0.9] )
		
		
		if plot:
			fig = plt.figure( figsize = (14,8) )
			ax = fig.add_axes( [0.05,0.15,0.45,0.8] , projection = "3d" )
			im_ld = ax.scatter( X[:,0] , X[:,1] , X[:,2] , s = 2. , c = ld , cmap = plt.cm.inferno , vmin = vmin_ld , vmax = vmax_ld )
			ax.set_title( "Local dimension (m = {})".format( ld.mean().round(2) ) )
			ax = fig.add_axes( [0.5,0.15,0.45,0.8] , projection = "3d" )
			im_th = ax.scatter( X[:,0] , X[:,1] , X[:,2] , s = 2. , c = th , cmap = plt.cm.inferno , vmin = vmin_th , vmax = vmax_th )
			ax.set_title( "Persistence" )
			
			cax = fig.add_axes( [ 0.12 , 0.1 , 0.3 , 0.05 ] )
			plt.colorbar( mappable = im_ld , cax = cax , orientation = "horizontal" )
			
			cax = fig.add_axes( [ 0.57 , 0.1 , 0.3 , 0.05 ] )
			plt.colorbar( mappable = im_th , cax = cax , orientation = "horizontal" )
			
			plt.show()
		print( "Test local dimension (Done)  " )
	except:
		print( "Test local dimension (Fail)  " )
##}}}


def test_mandelbrot( plot = True ):##{{{
	print( "Test Mandelbrot..." , end = "\r" )
	try:
		m = ckf.Mandelbrot( classic_set = "set0" , nx = 1000 , ny = 1000 )
		m.run()
		if plot:
			m.plot(True)
		print( "Test Mandelbrot. (Done)" )
	except:
		print( "Test Mandelbrot. (Fail)" )
##}}}

def test_julia( plot = True ):##{{{
	print( "Test Julia..." , end = "\r" )
	try:
		jul = ckf.Julia( classic_set = "set0" , nx = 1000 , ny = 1000 )
		jul.run()
		if plot:
			jul.plot(True)
		print( "Test Julia. (Done)" )
	except:
		print( "Test Julia. (Fail)" )
##}}}


def run_all_tests( plot = False ):##{{{
	test_lorenz63(plot)
	test_rossler(plot)
	test_lorenz84(plot)
	test_henon(plot)
	test_ikeda(plot)
	test_mira(plot)
	
	test_local_dimension(plot)
	
	test_mandelbrot(plot)
	test_julia(plot)
##}}}


#############
## Classes ##
#############

def figure_attractors():##{{{
	l63 = ck.Lorenz63()
	X0  = l63.orbit( np.linspace( 0 , 100 , 100 ) )[-1,:]
	X   = l63.orbit( np.linspace( 0 , 100 , 10000 ) , X0 = X0 )
	
	mir = ck.Mira()
	Y0  = mir.orbit( 10000 , X0 = np.array( [0.5,0.5] ) )
	mir = ck.Mira( size = 10000 )
	Y   = mir.orbit( 1000 , X0 = Y0 )[-1,:,:]
	
	
	fig = plt.figure( figsize = (14,7) )
	
	ax = fig.add_subplot( 1 , 2 , 1 , projection = "3d" )
	ax.plot( X[:,0] , X[:,1] , X[:,2] , color = "red" , linestyle = "-" )
	ax.set_xlabel( r"$x$" )
	ax.set_ylabel( r"$y$" )
	ax.set_zlabel( r"$z$" )
	ax.set_title( "Lorenz (1963) attractor" )
	
	ax = fig.add_subplot( 1 , 2 , 2 )
	ax.plot( Y[:,0] , Y[:,1] , color = "red" , linestyle = "" , marker = "." , markersize = 2 )
	ax.set_xlabel( r"$x$" )
	ax.set_ylabel( r"$y$" )
	ax.set_title( "Mira attractor" )
	
	fig.set_tight_layout(True)
	plt.savefig( "attractors.png" )
##}}}

def figure_snapshot_attractors():##{{{
	l63 = ck.Lorenz84( size = 5000 , F = "cyclic" )
	X0  = l63.orbit( np.linspace( 0 , 2 * 73 , 50 ) )[-1,:,:]
	X   = l63.orbit( np.linspace( 0 , 5 * 73 , 10000 ) , X0 = X0 )
	t_fall = int( 4    / 5 * 10000 )
	t_wint = int( 4.25 / 5 * 10000 )
	t_spri = int( 4.5  / 5 * 10000 )
	t_summ = int( 4.75 / 5 * 10000 )
	lX = []
	lX.append( X[t_fall,:,:] )
	lX.append( X[t_wint,:,:] )
	lX.append( X[t_spri,:,:] )
	lX.append( X[t_summ,:,:] )
	title = ["Fall","Winter","Spring","Summer"]
	color = ["orange" , "blue" , "green" , "red" ]
	
	fig = plt.figure( figsize = (12,12) )
	for i,X in enumerate(lX):
		ax = fig.add_subplot( 2 , 2 , i + 1 , projection = "3d" )
		ax.plot( X[:,0] , X[:,1] , X[:,2] , color = color[i] , linestyle = "" , marker = "." , markersize = 2 )
		ax.set_title( title[i] , fontsize = 20 )
	fig.set_tight_layout(True)
	plt.savefig("snapshot_attractors.png")
##}}}

def figure_local_dimension():##{{{
		l63 = ck.Lorenz63()
		X0  = l63.orbit( np.linspace( 0 , 100 , 100 ) )[-1,:]
		X   = l63.orbit( np.linspace( 0 , 100 , 10000 ) , X0 = X0 )
		
		ld,th = ck.localDimension( X , n_jobs = 6 , pareto_fit = "SDFC" )
		vmin_ld,vmax_ld = np.quantile( ld , [0.1,0.9] )
		vmin_th,vmax_th = np.quantile( th , [0.1,0.9] )
		
		
		fig = plt.figure( figsize = (14,8) )
		ax = fig.add_axes( [0.05,0.15,0.45,0.8] , projection = "3d" )
		im_ld = ax.scatter( X[:,0] , X[:,1] , X[:,2] , s = 2. , c = ld , cmap = plt.cm.inferno , vmin = vmin_ld , vmax = vmax_ld )
		ax.set_title( "Mean local dimension = {}".format( ld.mean().round(2) ) )
		
		ax = fig.add_axes( [0.5,0.15,0.45,0.8] , projection = "3d" )
		im_th = ax.scatter( X[:,0] , X[:,1] , X[:,2] , s = 2. , c = th , cmap = plt.cm.inferno , vmin = vmin_th , vmax = vmax_th )
		
		cax = fig.add_axes( [ 0.12 , 0.1 , 0.3 , 0.05 ] )
		cbar = plt.colorbar( mappable = im_ld , cax = cax , orientation = "horizontal" )
		cbar.set_label( "Local dimension" )
		
		cax = fig.add_axes( [ 0.57 , 0.1 , 0.3 , 0.05 ] )
		cbar = plt.colorbar( mappable = im_th , cax = cax , orientation = "horizontal" )
		cbar.set_label( "Persistence" )
		
		plt.savefig("local_dimension.png")
##}}}

def figure_mandelbrot_julia():
	
	m = ckf.Mandelbrot( classic_set = "set0" , nx = 2000 , ny = 2000 )
	m.run()
	jul = ckf.Julia( classic_set = "set2" , nx = 2000 , ny = 2000 )
	jul.run()
	
	fig = plt.figure( figsize = (14,7) )
	
	ax = fig.add_subplot( 1 , 2 , 1 )
	X,Y = np.meshgrid( m._x , m._y )
	ax.imshow( np.rot90(m.ratio) , cmap = plt.cm.hot_r )
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_title( "Mandelbrot set" )
	
	ax = fig.add_subplot( 1 , 2 , 2 )
	X,Y = np.meshgrid( jul._x , jul._y )
	ax.imshow( np.rot90(jul.ratio) , cmap = plt.cm.hot_r )
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_title( "Julia set" )
	
	fig.set_tight_layout(True)
	plt.savefig( "mandelbrot_julia.png" )

##########
## main ##
##########

if __name__ == "__main__":
	
	## Figure of attractors
	##=====================
	
#	figure_attractors()
#	figure_snapshot_attractors()
#	figure_local_dimension()
	figure_mandelbrot_julia()
	
##![Alt](/figures/)
##![Alt](/figures/)
	
	print("Done")


