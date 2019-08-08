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
	l63 = sk.Lorenz63()
	X0  = l63.orbit( np.linspace( 0 , 100 , 100 ) )[-1,:]
	X   = l63.orbit( np.linspace( 0 , 100 , 10000 ) , X0 = X0 )
	
	fig = plt.figure( figsize = (7,7) )
	ax = fig.add_subplot( 1 , 1 , 1 , projection = "3d" )
	ax.plot( X[:,0] , X[:,1] , X[:,2] , color = "red" , linestyle = "-" )
	fig.set_tight_layout(True)
	plt.show()
	
	
	print("Done")


