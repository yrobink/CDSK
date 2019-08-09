# -*- coding: utf-8 -*-

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

#mpl.rcParams['font.size'] = 30
#plt.rc('text',usetex=True)
#plt.rcParams['text.latex.unicode'] = True


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
			ax.set_title( "Mean local dimension = {}".format( ld.mean().round(2) ) )
			ax = fig.add_axes( [0.5,0.15,0.45,0.8] , projection = "3d" )
			im_th = ax.scatter( X[:,0] , X[:,1] , X[:,2] , s = 2. , c = th , cmap = plt.cm.inferno , vmin = vmin_th , vmax = vmax_th )
			
			cax = fig.add_axes( [ 0.12 , 0.1 , 0.3 , 0.05 ] )
			cbar = plt.colorbar( mappable = im_ld , cax = cax , orientation = "horizontal" )
			cbar.set_label( "Local dimension" )
			
			cax = fig.add_axes( [ 0.57 , 0.1 , 0.3 , 0.05 ] )
			cbar = plt.colorbar( mappable = im_th , cax = cax , orientation = "horizontal" )
			cbar.set_label( "Persistence" )
			
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



##########
## main ##
##########

if __name__ == "__main__":
	
	## Start by print version number
	print(ck.__version__)
	
	## Now tests
	run_all_tests()
	 
	
	print("Done")


