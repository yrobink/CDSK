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

def julia( c = complex(0.3,0.5) , delta = 2 , nxy = 1500 , maxit = 200 ):
	
	x = np.linspace( c.real - delta , c.real + delta , nxy )
	y = np.linspace( c.imag - delta , c.imag + delta , nxy )
	
	z = np.array( [ complex(x[i],y[j]) for i,j in itt.product(range(nxy),range(nxy)) ] )
	idx = np.zeros( (nxy**2) )
	
	for _ in range(maxit):
		z = z**2 + c
		idx[ np.abs(z) < 2 ] += 1
	
	
	res = idx.reshape( (nxy,nxy) ) / maxit
	
	cmap = plt.cm.hot_r
	mpl.image.imsave( "Julia.png" , np.rot90(res) , cmap = cmap )
#	fig = plt.figure( figsize = (10,10) )
#	ax = fig.add_subplot( 1 , 1 , 1 )
#	ax.imshow( np.rot90(res) , cmap = cmap , extent = [x.min(),x.max(),y.min(),y.max()] )
#	plt.tight_layout()
#	plt.show()


