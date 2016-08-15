# -*- coding: utf8 -*-
# # #
# Downscale PCMDI AR5 data to a pre-processed climatology
#  extent, resolution, reference system
#
# Author: Michael Lindgren (malindgren@alaska.edu)
# # #
import rasterio, os
# from pathos import multiprocessing
import numpy as np
import pandas as pd
from downscale import utils

class Baseline( object ):
	'''
	simple class to store the baseline arr rasters
	'''
	def __init__( self, filelist ):
		'''
		class for the baseline arr used as the template climatology
		to downscale the anomalies of the series to.

		Arguments:
		----------
		filelist = [list] of str paths to each of the 12 monthly climatology files.
				* must be in chronological order jan-dec.
		'''
		self.filelist = sorted( filelist )
		self.meta = rasterio.open( self.filelist[0] ).meta
		self.arrlist = [ rasterio.open( fn ).read( 1 ) for fn in self.filelist ]
		
		# make sure the data are sorted chronlogically
		# arrlist = self.arrlist
		# self.arrlist = self.sort_files()
		# self.rasters = [ rasterio.open(i).read(1) for i in self.arrlist ]

	# def sort_files( self ):
	# 	'''
	# 	ensure the files are in chronological order
	# 	'''
	# 	arrlist = self.arrlist
	# 	arrlist.sort()
	# 	self.arrlist = arrlist

class Dataset( object ):
	def __init__( self, fn, variable, model, scenario, project=None, units=None, metric=None, interp=False, ncpus=32, \
					method='linear', begin=None, end=None, *args, **kwargs):
		'''
		fn = [str] path to the xray supported dataset to be read in.
		variable = [str] abbreviation of variable name to extract from file
		model = [str] name of the model being read
		scenario = [str] name of the scenario being read
		units = [str] abbreviation of the units of the variable
		interp = [bool] if True interpolate across NA's using a spline. 
					if False (default) do nothing.interp=False, ncpus=32,
		ncpus = [ int ] number of cores to use if interp=True. default:2.
		'''
		import xarray as xr
		self.fn = fn
		self.ds = xr.open_dataset( self.fn )
		self.variable = variable
		self.model = model
		self.scenario = scenario
		self.begin = begin # year begin for self._open_dataset()
		self.end = end # year end
		
		if units != None:
			self.units = units
		else:
			self.units = 'units'
			
		if project != None:
			self.project = project
		else:
			self.project = 'project'

		if metric != None:
			self.metric = metric
		else:
			self.metric = 'metric'

		self.interp = interp
		self.ncpus = ncpus
		self.method = method
		self._rotated = False
		self._lonpc = None
		if interp:
			print( 'running interpolation across NAs' )
			_ = self.interp_na( )
	
	# def _open_dataset( self ):
	# 	ds = xr.open_dataset( self.fn )
	# 	if self.begin and self.end:
	# 		ds = ds.sel( time=slice( self.begin, self.end ) )
	# 	return ds
	def _calc_affine( self, *args, **kwargs ):
		import affine
		import numpy as np
		lat_shape, lon_shape = self.ds.dims[ 'lat' ], self.ds.dims[ 'lon' ]
		lonmin = self.ds.lon.min().data
		latmax = self.ds.lat.max().data

		# get the right upper left lon to the corner
		lon_arr = np.array([-180.0, 0.0])
		idx = (np.abs(lon_arr - lonmin)).argmin()
		lonmin = lon_arr[ idx ]
		
		# HARDWIRED FOR GLOBAL DATA
		latmax = 90.0

		# get the resolution in both directions
		lat_res = 180.0 / lat_shape
		lon_res = 360.0 / lon_shape
		return affine.Affine( lon_res, 0.0, lonmin, 0.0, -lat_res, latmax )
	# def _calc_affine( self, *args, **kwargs ):
	# 	# POTENTIALLY REMOVE THIS FROM HERE?  IDK WHERE THIS IS BEST FIT
	# 	'''
	# 	this assumes 0-360 longitude-ordering (pacific-centered)
	# 	and WGS84 LatLong (Decimal Degrees). EPSG:4326.
	# 	'''
	# 	import affine
	# 	lat_shape, lon_shape = self.ds.dims[ 'lat' ], self.ds.dims[ 'lon' ]
	# 	lonmin = self.ds.lon.min().data
	# 	latmax = self.ds.lat.max().data
	# 	lat_res = 180.0 / lat_shape
	# 	lon_res = 360.0 / lon_shape
	# 	return affine.Affine( lon_res, 0.0, lonmin, 0.0, -lat_res, latmax )
	@staticmethod
	def rotate( dat, lons, to_pacific=False ):
		'''rotate longitudes in WGS84 Global Extent'''
		if to_pacific == True:
			# to 0 - 360
			dat, lons = utils.shiftgrid( 0., dat, lons )
		elif to_pacific == False:
			# to -180.0 - 180.0 
			dat, lons = utils.shiftgrid( 180., dat, lons, start=False )
		else:
			raise AttributeError( 'to_pacific must be boolean True:False' )
		return dat, lons
	@staticmethod
	def wrap( d ):
		return utils.xyz_to_grid( **d )
	def interp_na( self ):
		'''
		np.float32
		method = [str] one of 'cubic', 'near', 'linear'

		return a list of dicts to pass to the xyz_to_grid hopefully in parallel
		'''
		from copy import copy
		# from pathos import multiprocessing
		# from multiprocessing import Pool
		# from pathos.mp_map import mp_map
		import pandas as pd
		import numpy as np
		# remove the darn scientific notation
		np.set_printoptions( suppress=True )
		output_dtype = np.float32
		
		# if 0-360 leave it alone
		if ( self.ds.lon > 200.0 ).any() == True:
			dat, lons = self.ds[ self.variable ].data, self.ds.lon
			self._lonpc = lons
		else:
			# greenwich-centered rotate to 0-360 for interpolation across pacific
			dat, lons = self.rotate( self.ds[ self.variable ].values, self.ds.lon, to_pacific=True )
			self._rotated = True # update the rotated attribute
			self._lonpc = lons

		# mesh the lons and lats and unravel them to 1-D
		xi, yi = np.meshgrid( self._lonpc, self.ds.lat.data )
		lo, la = [ i.ravel() for i in (xi,yi) ]

		# setup args for multiprocessing
		df_list = [ pd.DataFrame({ 'x':lo, 'y':la, 'z':d.ravel() }).dropna( axis=0, how='any' ) for d in dat ]

		args = [ {'x':np.array(df['x']), 'y':np.array(df['y']), 'z':np.array(df['z']), \
				'grid':(xi,yi), 'method':self.method, 'output_dtype':output_dtype } for df in df_list ]

		# def wrap( d ):
		# 	return 

		print( 'processing cru re-gridding in serial due to multiprocessing issues...' )
		dat = np.array([ utils.xyz_to_grid( **i ) for i in args ])
		# pool = multiprocessing.Pool( self.ncpus )
		# dat = np.array( pool.map( self.wrap, args ) )
		# pool.close()
		# pool.join()

		lons = self._lonpc
		if self._rotated == True: # rotate it back
			dat, lons = self.rotate( dat, lons, to_pacific=False )
				
		# place back into a new xarray.Dataset object for further processing
		self.ds = self.ds.update( { self.variable:( ['time','lat','lon'], dat ) } )
		print( 'ds interpolated updated into self.ds' )
		return 1