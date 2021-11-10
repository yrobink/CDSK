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
import matplotlib.pyplot as plt
import matplotlib.gridspec as mplg

import CDSK as ck
import CDSK.fractal as ckf


####################
## mpl parameters ##
####################


mpl.rcdefaults()
mpl.rcParams['font.size'] = 7
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['lines.markersize'] = 1
mpl.rcParams['patch.linewidth'] = 0.5
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['ytick.major.width'] = 0.5


###############
## Functions ##
###############

def test_lorenz63():##{{{
	print( "Test Lorenz63..." , end = "\r" )
	try:
		l63 = ck.Lorenz63()
		X0  = l63.orbit( np.linspace( 0 , 100 , 100 ) )[-1,:]
		X   = l63.orbit( np.linspace( 0 , 100 , 10000 ) , X0 = X0 )
		
		mm  = 1. / 25.4
		fig = plt.figure( figsize = (90*mm,90*mm) )
		ax = fig.add_subplot( 1 , 1 , 1 , projection = "3d" )
		ax.plot( X[:,0] , X[:,1] , X[:,2] , color = "red" , linestyle = "-" )
		plt.subplots_adjust( left = 0 , right = 0.95 , bottom = 0 , top = 1 )
		plt.savefig( "output/test_lorenz63.png" , dpi = 600 )
		print( "Test Lorenz63 (Done)  " )
	except:
		print( "Test Lorenz63 (Fail)  " )
##}}}

def test_rossler():##{{{
	print( "Test Rossler..." , end = "\r" )
	try:
		ross = ck.Rossler()
		X0  = ross.orbit( np.linspace( 0 , 100 , 100 ) )[-1,:]
		X   = ross.orbit( np.linspace( 0 , 300 , 10000 ) , X0 = X0 )
		
		mm  = 1. / 25.4
		fig = plt.figure( figsize = (90*mm,90*mm) )
		ax = fig.add_subplot( 1 , 1 , 1 , projection = "3d" )
		ax.plot( X[:,0] , X[:,1] , X[:,2] , color = "red" , linestyle = "-" )
		plt.subplots_adjust( left = 0 , right = 0.95 , bottom = 0 , top = 1 )
		plt.savefig( "output/test_rossler.png" , dpi = 600 )
		print( "Test Rossler (Done)  " )
	except:
		print( "Test Rossler (Fail)  " )
##}}}

def test_lorenz84():##{{{
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
		
		mm  = 1. / 25.4
		fig = plt.figure( figsize = (90*mm,90*mm) )
		for i,X in enumerate(lX):
			ax = fig.add_subplot( 2 , 2 , i + 1 , projection = "3d" )
			ax.plot( X[:,0] , X[:,1] , X[:,2] , color = "red" , linestyle = "" , marker = "." )
			ax.set_title( title[i] )
		plt.subplots_adjust( left = 0 , right = 0.92 , bottom = 0.05 , top = 0.95 )
		plt.savefig( "output/test_lorenz84.png" , dpi = 600 )
		print( "Test Lorenz84 (Done)  " )
	except:
		print( "Test Lorenz84 (Fail)  " )
##}}}

def test_henon():##{{{
	print( "Test Henon..." , end = "\r" )
	try:
		henon = ck.Henon()
		X   = henon.orbit( 10000 , X0 = np.array( [0.5,0.5] ) )
		
		mm  = 1. / 25.4
		fig = plt.figure( figsize = (90*mm,90*mm) )
		ax = fig.add_subplot( 1 , 1 , 1 )
		ax.plot( X[:,0] , X[:,1] , color = "red" , linestyle = "" , marker = "." )
		plt.tight_layout()
		plt.savefig( "output/test_henon.png" , dpi = 600 )
		print( "Test Henon (Done)  " )
	except:
		print( "Test Henon (Fail)  " )
##}}}

def test_ikeda():##{{{
	print( "Test Ikeda..." , end = "\r" )
	try:
		ike = ck.Ikeda()
		X   = ike.orbit( 10000 , X0 = np.array( [0.5,0.5] ) )
		
		mm  = 1. / 25.4
		fig = plt.figure( figsize = (90*mm,90*mm) )
		ax = fig.add_subplot( 1 , 1 , 1 )
		ax.plot( X[:,0] , X[:,1] , color = "red" , linestyle = "" , marker = "." )
		plt.tight_layout()
		plt.savefig( "output/test_ikeda.png" , dpi = 600 )
		print( "Test Ikeda (Done)  " )
	except:
		print( "Test Ikeda (Fail)  " )
##}}}

def test_mira():##{{{
	print( "Test Mira..." , end = "\r" )
	try:
		mir = ck.Mira()
		X0  = mir.orbit( 10000 , X0 = np.array( [0.5,0.5] ) )
		mir = ck.Mira( size = 10000 )
		X   = mir.orbit( 1000 , X0 = X0 )[-1,:,:]
		
		mm  = 1. / 25.4
		fig = plt.figure( figsize = (90*mm,90*mm) )
		ax = fig.add_subplot( 1 , 1 , 1 )
		ax.plot( X[:,0] , X[:,1] , color = "red" , linestyle = "" , marker = "." )
		plt.tight_layout()
		plt.savefig( "output/test_mira.png" , dpi = 600 )
		print( "Test Mira (Done)  " )
	except:
		print( "Test Mira (Fail)  " )
##}}}


def test_local_dimension():##{{{
	print( "Test local dimension..." , end = "\r" )
	try:
		## Dyn. indexes at points defined by X
		l63 = ck.Lorenz63( size = 1000 )
		X01 = l63.orbit( np.linspace( 0 , 100 , 10000 ) )[-1,:,:]
		
		## Two reference trajectories
		l63 = ck.Lorenz63()
		y0  = l63.orbit( np.linspace( 0 , 100 , 100 ) )[-1,:]
		y1  = y0 + np.random.normal( scale = 0.1 , size = 3 )
		Y0  = l63.orbit( np.linspace( 0 , 100 , 20000 ) , X0 = y0 )
		Y1  = l63.orbit( np.linspace( 0 , 100 , 20000 ) , X0 = y1 )
		
		## Merge
		X = np.stack( (X01,X01) , -1 )
		Y = np.stack( (Y0,Y1) , -1 )
		
		ld,theta,alpha,shape,where = ck.dynamical_local_indexes( X , Y = Y , ld_fit = "mean" , n_jobs = 6 )
		
		fig  = plt.figure()
		grid = mplg.GridSpec( 6 , 3 )
		
		## Local dim
		ii = [0,1,0]
		jj = [0,1,1]
		for i in range(3):
			ax = fig.add_subplot( grid[i,0] , projection = "3d" )
			im = ax.scatter( X01[:,0] , X01[:,1] , X01[:,2] , c = ld[:,ii[i],jj[i]] , cmap = plt.cm.inferno , vmin = 1 , vmax = 3 )
		
		gax = grid[4,0].subgridspec( 1 , 3 , width_ratios = [0.05,1,0.05] )
		cax = fig.add_subplot( gax[0,1] )
		plt.colorbar( mappable = im , cax = cax , orientation = "horizontal" , label = "Local dimension" )
		
		## Theta
		for i in range(3):
			ax = fig.add_subplot( grid[i,2] , projection = "3d" )
			im = ax.scatter( X01[:,0] , X01[:,1] , X01[:,2] , c = theta[:,ii[i],jj[i]] , cmap = plt.cm.inferno , vmin = 0 , vmax = 0.5 )
		
		gax = grid[4,2].subgridspec( 1 , 3 , width_ratios = [0.05,1,0.05] )
		cax = fig.add_subplot( gax[0,1] )
		plt.colorbar( mappable = im , cax = cax , orientation = "horizontal" , label = r"$\theta$" )
		
		
		## Figsize
		mm    = 1. / 25.4
		width = 110*mm
		width_cax = width / 2
		widths = [width_cax,10*mm,width_cax]
		
		height_cax = width_cax * ax.get_data_ratio() * 0.77
		heights = [height_cax for _ in range(3)] + [3*mm] + [5*mm] + [10*mm]
		height = np.sum(heights)
		
		grid.set_height_ratios(heights)
		grid.set_width_ratios(widths)
		fig.set_figheight(height)
		fig.set_figwidth(width)
		
		plt.subplots_adjust( left = 0 , bottom = 0 , right = 1 , top = 1 , wspace = 0 , hspace = 0 )
		plt.savefig( "output/dyn_index.png" , dpi = 600 )
		plt.close(fig)
		print( "Test local dimension (Done)  " )
	except:
		print( "Test local dimension (Fail)  " )
##}}}


def test_mandelbrot():##{{{
	print( "Test Mandelbrot..." , end = "\r" )
	try:
		m = ckf.Mandelbrot( classic_set = "set0" , nx = 1000 , ny = 1000 )
		m.run()
		
		##
		X,Y = np.meshgrid( m.x , m.y )
		Z   = m.ratio
		
		fig = plt.figure()
		ax  = fig.add_subplot(1,1,1)
		ax.pcolormesh( X , Y , Z.T , cmap = plt.cm.hot_r , shading = "nearest" )
		
		mm    = 1. / 25.4
		width = 100*mm
		height = width * ax.get_data_ratio()
		
		fig.set_figwidth(width)
		fig.set_figheight(height)
		
		plt.subplots_adjust( left = 0 , bottom = 0 , right = 1 , top = 1 )
		plt.savefig( "output/test_mandelbrot.png" , dpi = 600 )
		
		print( "Test Mandelbrot. (Done)" )
	except:
		print( "Test Mandelbrot. (Fail)" )
##}}}

def test_julia():##{{{
	print( "Test Julia..." , end = "\r" )
	try:
		jul = ckf.Julia( classic_set = "set0" , nx = 1000 , ny = 1000 )
		jul.run()
		
		##
		X,Y = np.meshgrid( jul.x , jul.y )
		Z   = jul.ratio
		
		fig = plt.figure()
		ax  = fig.add_subplot(1,1,1)
		ax.pcolormesh( X , Y , Z.T , cmap = plt.cm.hot_r , shading = "nearest" )
		
		mm    = 1. / 25.4
		width = 100*mm
		height = width * ax.get_data_ratio()
		
		fig.set_figwidth(width)
		fig.set_figheight(height)
		
		plt.subplots_adjust( left = 0 , bottom = 0 , right = 1 , top = 1 )
		plt.savefig( "output/test_julia.png" , dpi = 600 )
		print( "Test Julia. (Done)" )
	except:
		print( "Test Julia. (Fail)" )
##}}}


def run_all_tests():##{{{
	if not os.path.isdir("output"):
		os.makedirs("output")
	test_lorenz63()
	test_rossler()
	test_lorenz84()
	test_henon()
	test_ikeda()
	test_mira()
	
	test_local_dimension()
	
	test_mandelbrot()
	test_julia()
##}}}


#############
## Classes ##
#############



##########
## main ##
##########

if __name__ == "__main__":
	
	## Start by print version number
	print(ck.__version__)
	
	## Now tests
	run_all_tests()
	
	
	print("Done")


