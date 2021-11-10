# -*- coding: utf-8 -*-

#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

###############
## Libraries ##
##{{{

import sys,os
import itertools as itt

## Scientific libraries
##=====================
import numpy as np
import xarray as xr
import CDSK as cs


## Plot libraries ##
##==================
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as mplgrid

##}}}
###############

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
	
	np.random.seed(42)
	
	f_data = "data_l63.nc"
	if os.path.exists(f_data):
		df = xr.open_dataset(f_data)
	else:
		## Data
		l63 = cs.Lorenz63()
		t0  = np.arange( 0 , 10 , 0.01 )
		t   = np.arange( 0 , 100 , 0.01 )
		Y0 = l63.orbit( t0 )[-1,:]
		Y  = l63.orbit( t , Y0 )
		
		n_snap = 1000
		l63 = cs.Lorenz63( size = n_snap )
		X = l63.orbit( t )[-1,:,:]
		
		## Save
		R = ["R{:{fill}{align}{n}}".format(i,fill=0,align=">",n=3) for i in range(n_snap)]
		X     = xr.DataArray( X     , dims = ["ref","axes"]  , coords = [R,["x","y","z"]] )
		Y     = xr.DataArray( Y     , dims = ["time","axes"] , coords = [t,["x","y","z"]] )
		df    = xr.Dataset( { "X" : X , "Y" : Y } )
		df.to_netcdf(f_data)
	
	## Arguments
	kwargs = { "ld_fit" : "SDFC" }
	
	## Pairwise projections
	X2 = np.zeros( (df.ref.size,2,3) )
	X2[:,:,0] = df.X[:,[0,1]]
	X2[:,:,1] = df.X[:,[0,2]]
	X2[:,:,2] = df.X[:,[1,2]]
	Y2 = np.zeros( (df.time.size,2,3) )
	Y2[:,:,0] = df.Y[:,[0,1]]
	Y2[:,:,1] = df.Y[:,[0,2]]
	Y2[:,:,2] = df.Y[:,[1,2]]
	ld2,theta2,alpha2,shp2,where2 = cs.dynamical_local_indexes( X2 , Y2 , **kwargs )
	
	## 3d attractor
	X3 = df.X.values
	Y3 = df.Y.values
	ld3,theta3,alpha3,shp3,where3 = cs.dynamical_local_indexes( X3 , Y3 , **kwargs )
	
	## Figures
	l_cmap = [plt.cm.inferno,plt.cm.turbo,plt.cm.binary,plt.cm.coolwarm]
	nrow,ncol = 4,3
	grid = mplgrid.GridSpec( nrow , ncol , height_ratios = [1,1,1,0.05] )
	
	for c,c3d,title,cmap in zip([ld2,theta2,alpha2,shp2],[ld3,theta3,alpha3,shp3],["Local Dimension","Theta","Alpha","Shape"],l_cmap):
		fig = plt.figure( figsize = (12,12) )
		
		vmin2,vmax2 = np.nanquantile( c   , [0.1,0.9] )
		vmin3,vmax3 = np.nanquantile( c3d , [0.1,0.9] )
		vmin = min(vmin2,vmin3)
		vmax = max(vmax2,vmax3)
		if title == "Shape":
			vminmax = max(abs(vmin),abs(vmax))
			vmin = -vminmax
			vmax = vminmax
		
		kwargs = { "vmin" : vmin , "vmax" : vmax , "cmap" : cmap , "s" : 5 }
		for i,j in itt.combinations_with_replacement(range(3),2):
			ax = fig.add_subplot( grid[i,j] )
			
			if i == j:
				im = ax.scatter( X2[:,0,i] , X2[:,1,i] , c = c[:,i,j] , **kwargs )
			else:
				ax.plot( c[:,i,j] , color = "blue" )
			
			ax.set_title( "Mean = {:.2f}".format( np.mean(c[:,i,j]) ) )
		
		cax = fig.add_subplot( grid[3,:] )
		plt.colorbar( mappable = im , cax = cax , orientation = "horizontal" , label = title )
		
		ax  = fig.add_subplot( grid[2,0] , projection = "3d" )
		ax.scatter( X3[:,0] , X3[:,1] , X3[:,2] , c = c3d[:,0,0] , **kwargs )
		ax.set_title( "Mean = {:.2f}".format( np.mean(c3d) ) )
		
		plt.tight_layout()
		plt.savefig( "fig_{}.png".format(title) )
		plt.close(fig)
	
	
	
	print("Done")

