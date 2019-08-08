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


import CDSK as sk


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
		l63 = sk.Lorenz63()
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
		ross = sk.Rossler()
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
		l63 = sk.Lorenz84( size = 1000 , F = "cyclic" )
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
		henon = sk.Henon()
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
		ike = sk.Ikeda()
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
		mir = sk.Mira()
		X0  = mir.orbit( 10000 , X0 = np.array( [0.5,0.5] ) )
		mir = sk.Mira( size = 10000 )
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


def run_all_tests( plot = False ):##{{{
	test_lorenz63(plot)
	test_rossler(plot)
	test_lorenz84(plot)
	test_henon(plot)
	test_ikeda(plot)
	test_mira(plot)
##}}}


#############
## Classes ##
#############



##########
## main ##
##########

if __name__ == "__main__":
	
	## Start by print version number
	print(sk.__version__)
	
	## Now tests
	run_all_tests()
	
	print("Done")


