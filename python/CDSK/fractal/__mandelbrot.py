# -*- coding: utf-8 -*-

###############
## Libraries ##
###############

import numpy as np
import itertools as itt
import matplotlib as mpl


###############
## Functions ##
###############

def mandelbrot( box = "set0" , nx = 1500 , ny = 1500 , maxit = 200 ):
	
	if type(box) == tuple:
		x,y = box
	else:
		if box == "set1":
			x = np.linspace( -0.62 , -0.42 , nx )
			y = np.linspace( -0.7 , -0.5 , ny )
		elif box == "set2":
			x = np.linspace( -0.60 , -0.595 , nx )
			y = np.linspace( -0.665 , -0.66 , ny )
		elif box == "set3":
			x = np.linspace( -0.57 , -0.56 , nx )
			y = np.linspace( -0.646 , -0.636 , ny )
		elif box == "set4":
			x = np.linspace( -0.5663 , -0.5656 , nx )
			y = np.linspace( -0.6394 , -0.6387 , ny )
		else:
			x = np.linspace( -1.5 , -0.5 , nx )
			y = np.linspace( -0.5 , 0.5 , ny )
	
	c = np.array( [ complex(x[i],y[j]) for i,j in itt.product(range(nx),range(ny)) ] )
	z = np.zeros( (nx*ny) )
	idx = np.zeros( (nx*ny) )
	
	for _ in range(maxit):
		z = z**2 + c
		idx[ np.abs(z) < 2 ] += 1
	
	
	res = idx.reshape( (nx,ny) ) / maxit
	
	fig = plt.figure( figsize = (nx / ny * 10 , 10 ) )
	ax = fig.add_subplot( 1 , 1 , 1 )
	ax.imshow( np.rot90(res) , cmap = plt.cm.hot_r , extent = [x.min(),x.max(),y.min(),y.max()] )
	ax.set_axis_off()
	plt.tight_layout()
	plt.show()
#	plt.savefig( "Test.png" )




