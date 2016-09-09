# -*- coding: utf8 -*-
# # #
# Downscale PCMDI AR5 data to a pre-processed climatology
#  extent, resolution, reference system
#
# Author: Michael Lindgren (malindgren@alaska.edu)
# # #
import rasterio, os
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
		self.arrlist = ( rasterio.open( fn ).read( 1 ) for fn in self.filelist )

# this may be a ridiculous class considering the flexibility of its big brother.
class BaselineTS( object ):
	'''
	class to read / sort / store the mean downscaled data series
	used in min/max delta-style downscaling.  This is different to preserve
	the relationships between the max/mean/min variables that diverge using
	the standard DeltaDownscale approach. This will be a full timeseries of
	pre-downscaled `tas` data...  
	'''
	def __init__( self, filelist )


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
		northup = [bool] if True, flip the earth using np.flipud if False leave it alone
		'''
		import xarray as xr
		self.fn = fn
		self.ds = xr.open_dataset( self.fn )
		self.variable = variable
		self.model = model
		self.scenario = scenario
		self.begin = begin # year begin
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

		# slice to the years we want if given
		if self.begin != None and self.end != None:
			self.ds = self.ds.sel( time=slice( str( self.begin ), str( self.end ) ) )

		self.interp = interp
		self.ncpus = ncpus
		self.method = method
		self._rotated = False
		self._lonpc = None
		if interp:
			print( 'running interpolation across NAs' )
			_ = self.interp_na( )
	
	@staticmethod
	def transform_from_latlon( lat, lon ):
		''' simple way to make an affine transform from lats and lons coords '''
		from affine import Affine
		lat = np.asarray( lat )
		lon = np.asarray( lon )
		trans = Affine.translation(lon[0], lat[0])
		scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
		return trans * scale
	# @staticmethod
	# def transform_from_latlon( lat, lon ):
	# 	''' simple way to make an affine transform from lats and lons coords '''
	# 	from affine import Affine
	# 	lat = np.asarray( lat )
	# 	lon = np.asarray( lon )
	# 	if (np.max( lat ) - 90) < np.abs( np.mean( np.diff( lat ) ) ):
	# 		lat_max = 90.0
	# 	else:
	# 		lat_max = np.max( lat )

	# 	# set the lonmax to the corner.
	# 	lon_arr = np.array([-180.0, 0.0 ])
	# 	idx = (np.abs(lon_arr - np.min( lon ) ) ).argmin()
	# 	lon_max = lon_arr[ idx ]

	# 	trans = Affine.translation(lon_max, lat_max)
	# 	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
	# 	return trans * scale
	def _calc_affine( self ):
		return self.transform_from_latlon( self.ds.lat, self.ds.lon )
	def _northup( self, latitude='lat' ):
		''' this works only for global grids to be downscaled flips it northup '''
		if self.ds[ latitude ][0].data < 0: # meaning that south is north globally
			self.ds[ latitude ] = np.flipud( self.ds[ latitude ] )
			# flip each slice of the array and make a new one
			flipped = np.array( [ np.flipud( arr ) for arr in self.ds[ self.variable ].data ] )
			self.ds[ self.variable ] = (('time', 'lat', 'lon' ), flipped )
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
		import pandas as pd
		import numpy as np
		from pathos.mp_map import mp_map

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
		
		# # # # USE MLAB's griddata which we _can_ parallelize
		def wrap( d ):
			''' simple wrapper around utils.xyz_to_grid for mp_map'''
			x = np.array( d['x'] )
			y = np.array( d['y'] )
			z = np.array( d['z'] )
			xi, yi = d['grid']
			return utils.xyz_to_grid( x, y, z, (xi,yi), interp='linear' )
		# # # # 

		try:
			print( 'processing interpolation to convex hull in parallel using {} cpus.'.format( self.ncpus ) )
			dat_list = mp_map( wrap, args, nproc=self.ncpus )
			dat_list = [ i.data for i in dat_list ] # drop the output mask
			dat = np.array( dat_list )
		except:
			print( 'processing cru re-gridding in serial due to multiprocessing issues...' )
			dat = np.array([ wrap( **i ) for i in args ])

		lons = self._lonpc
		if self._rotated == True: # rotate it back
			dat, lons = self.rotate( dat, lons, to_pacific=False )
				
		# place back into a new xarray.Dataset object for further processing
		self.ds = self.ds.update( { self.variable:( ['time','lat','lon'], dat ) } )
		print( 'ds interpolated updated into self.ds' )
		return 1