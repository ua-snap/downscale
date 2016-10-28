# -*- coding: utf8 -*-
# # #
# Downscale PCMDI AR5 data to a pre-processed climatology
#  extent, resolution, reference system
#
# Author: Michael Lindgren (malindgren@alaska.edu)
# # #

import rasterio, os, copy
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from downscale import utils

class DeltaDownscale( object ):
	def __init__( self, baseline, clim_begin, clim_end, historical, future=None,
				downscaling_operation='add', level=None, level_name=None, 
				mask=None, mask_value=0,ncpus=32, src_crs={'init':'epsg:4326'}, 
				src_nodata=-9999.0, dst_nodata=None, post_downscale_function=None, 
				varname=None, modelname=None, anom=False, resample_type='bilinear', 
				fix_clim=False, interp=False, find_bounds=False, aoi_mask=None, *args, **kwargs ):
		
		'''
		simple delta downscaling
		
		Arguments:

		----------
		baseline = []
		clim_begin = []
		clim_end = []
		historical = []
		future = []
		level = []
		level_name = []
		
		Returns:
		--------
		
		'''
		self.historical = historical
		self.future = future
		self.baseline = baseline
		self.clim_begin = clim_begin
		self.clim_end = clim_end
		self.downscaling_operation = downscaling_operation
		self.level = level
		self.level_name = level_name
		self.mask = mask
		self.mask_value = mask_value
		self.ncpus = ncpus
		self.varname = varname
		self.modelname = modelname
		self.anom = anom
		self.resample_type = resample_type
		self.affine = self.historical._calc_affine()
		self.src_crs = src_crs
		self.src_nodata = src_nodata
		self.dst_nodata = dst_nodata
		self.post_downscale_function = post_downscale_function
		self.fix_clim = fix_clim
		self.interp = interp
		self.find_bounds = find_bounds
		self.aoi_mask = aoi_mask
		self.utils = utils

		# interpolate across space GCLL/PCLL args
		self._rotated = False
		self._lonpc = None

		# empty attributes to calculate
		self.anomalies = None
		self.climatology = None
		self.ds = None
		self._concat_nc()

		# fix pr climatologies if desired
		if fix_clim == True:
			print( 'fixing high/low values -- {}...'.format( self.varname ) )
			self.interp = True # force True 
			
			if self.aoi_mask is not None: # hairy
				mask = self.aoi_mask.mask
			else:
				mask = None

			self._calc_climatolgy()
			self._fix_clim( aoi_mask=mask, find_bounds=self.find_bounds )
			
			# interpolate clims across space
			self._interp_na_fix_clim()
			
			# fix the ds values -- will be interped below...
			self._fix_ds( aoi_mask=mask, find_bounds=self.find_bounds )

		if self.interp == True:
			print( 'running interpolation across NAs -- base resolution' )
			_ = self.interp_na( )

		# calculate climatology if fix_clim == False
		if self.fix_clim == False:
			self._calc_climatolgy()
			
		# calculate anomalies with the new climatology values
		self._calc_anomalies()

	def _concat_nc( self ):
		if self.historical and self.future:
			ds = xr.concat([ self.historical.ds, self.future.ds ], dim='time' )
		else:
			ds = self.historical.ds
		ds = ds[ self.historical.variable ]
		if self.level:
			# levidx, = np.where( ds[ self.level_name ] == self.level )
			# ds = ds[ :, levidx[0], ... ]
			ds = ds[ :, self.level, ... ]
		self.ds = ds
	def _calc_climatolgy( self ):
		'''slice / aggregate to climatology using mean'''
		try:
			climatology = self.ds.sel( time=slice( self.clim_begin, self.clim_end ) )
			self.climatology = climatology.groupby( 'time.month' ).mean( 'time' )
		except Exception:
			raise AttributeError( 'non-overlapping climatology period and series' )
	def _calc_anomalies( self ):
		''' calculate simple absolute or relative anomalies depending on variable '''
		if self.downscaling_operation == 'add':
			anomalies = self.ds.groupby( 'time.month' ) - self.climatology
		elif self.downscaling_operation == 'mult':
			anomalies = self.ds.groupby( 'time.month' ) / self.climatology
		else:
			NameError( '_calc_anomalies (ar5): value of downscaling_operation must be "add" or "mult" ' )

		# slice back to times we want
		if self.historical != None and self.future != None:
			self.anomalies = anomalies.sel( time=self.future.ds.time )
		else:
			self.anomalies = anomalies.sel( time=self.historical.ds.time )
	def _fix_clim( self, aoi_mask, find_bounds=False ):
		''' fix values in precip data '''
		if find_bounds == True:
			bound_mask = find_boundary( self.climatology[ 0, ... ].data )
			for idx in range( self.climatology.shape[0] ):
				arr = self.climatology[ idx, ... ].data
				arr = correct_boundary( arr, bound_mask, aoi_mask=aoi_mask )
				self.climatology[ idx, ... ].data = correct_inner( arr, bound_mask, aoi_mask=aoi_mask )

		elif find_bounds == False:
			for idx in range( self.climatology.shape[0] ):
				arr = self.climatology[ idx, ... ].data
				self.climatology[ idx, ... ].data = correct_values( arr, aoi_mask=aoi_mask )
		else:
			ValueError( 'find_bounds arg is boolean only' )
	def _fix_ds( self, aoi_mask, find_bounds=False ):
		''' fix high/low values in precip data '''
		if find_bounds == True:
			bound_mask = find_boundary( self.ds[ 0, ... ].data )
			for idx in range( self.ds.shape[0] ):
				arr = self.ds[ idx, ... ].data
				arr = correct_boundary( arr, bound_mask, aoi_mask=self.aoi_mask.mask )
				self.ds[ idx, ... ].data = correct_inner( arr, bound_mask, aoi_mask=aoi_mask )

		elif find_bounds == False:
			for idx in range( self.ds.shape[0] ):
				arr = self.ds[ idx, ... ].data
				self.ds[ idx, ... ].data = correct_values( arr, aoi_mask=aoi_mask )
		else:
			ValueError( 'find_bounds arg is boolean only' )
	@staticmethod
	def wrap( d ):
		''' simple wrapper around utils.xyz_to_grid for mp_map '''
		x = np.array( d['x'] )
		y = np.array( d['y'] )
		z = np.array( d['z'] )
		xi, yi = d['grid']
		return utils.xyz_to_grid( x, y, z, (xi,yi), interp='linear' )
	def interp_na( self ):
		'''
		np.float32
		method = [str] one of 'cubic', 'near', 'linear'

		return a list of dicts to pass to the xyz_to_grid in parallel
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
			dat, lons = self.ds.data, self.ds.lon
			self._lonpc = lons
		else:
			# greenwich-centered rotate to 0-360 for interpolation across pacific
			dat, lons = self.utils.rotate( self.ds.values, self.ds.lon, to_pacific=True )
			self._rotated = True # update the rotated attribute
			self._lonpc = lons

		# mesh the lons and lats and unravel them to 1-D
		xi, yi = np.meshgrid( self._lonpc, self.ds.lat.data )
		lo, la = [ i.ravel() for i in (xi,yi) ]

		# setup args for multiprocessing
		df_list = [ pd.DataFrame({ 'x':lo, 'y':la, 'z':d.ravel() }).dropna( axis=0, how='any' ) for d in dat ]

		args = [ {'x':np.array(df['x']), 'y':np.array(df['y']), 'z':np.array(df['z']), \
				'grid':(xi,yi), 'method':self.historical.method, 'output_dtype':output_dtype } for df in df_list ]
		
		print( 'processing interpolation to convex hull in parallel using {} cpus.'.format( self.ncpus ) )
		dat_list = mp_map( self.wrap, args, nproc=self.ncpus )
		dat_list = [ np.array(i) for i in dat_list ] # drop the output mask
		dat = np.array( dat_list )

		lons = self._lonpc
		if self._rotated == True: # rotate it back
			dat, lons = self.utils.rotate( dat, lons, to_pacific=False )
				
		# place back into a new xarray.Dataset object for further processing
		# self.ds = self.ds.update( { self.historical.variable:( ['time','lat','lon'], dat ) } )
		self.ds.data = dat
		print( 'ds interpolated updated into self.ds' )
		return 1
	def _interp_na_fix_clim( self ):
		'''
		np.float32
		method = [str] one of 'cubic', 'near', 'linear'

		return a list of dicts to pass to the xyz_to_grid in parallel
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
			dat, lons = self.climatology.data, self.ds.lon
			self._lonpc = lons
			self._rotated = False
		else:
			# greenwich-centered rotate to 0-360 for interpolation across pacific
			dat, lons = self.utils.rotate( self.climatology.values, self.ds.lon, to_pacific=True )
			self._rotated = True # update the rotated attribute
			self._lonpc = lons

		# mesh the lons and lats and unravel them to 1-D
		xi, yi = np.meshgrid( self._lonpc, self.ds.lat.data )
		lo, la = [ i.ravel() for i in (xi,yi) ]

		# setup args for multiprocessing
		df_list = [ pd.DataFrame({ 'x':lo, 'y':la, 'z':d.ravel() }).dropna( axis=0, how='any' ) for d in dat ]

		args = [ {'x':np.array(df['x']), 'y':np.array(df['y']), 'z':np.array(df['z']), \
				'grid':(xi,yi), 'method':self.historical.method, 'output_dtype':output_dtype } for df in df_list ]
		
		print( 'processing interpolation to convex hull in parallel using {} cpus. -- CLIMATOLOGY'.format( self.ncpus ) )
		dat_list = mp_map( self.wrap, args, nproc=self.ncpus )
		dat_list = [ np.array(i) for i in dat_list ] # drop the output mask
		dat = np.array( dat_list )

		lons = self._lonpc
		if self._rotated == True: # rotate it back
			dat, lons = self.utils.rotate( dat, lons, to_pacific=False )
			self._rotated = False # reset it now that its back
				
		# place back into a new xarray.Dataset object for further processing
		# self.ds = self.ds.update( { self.historical.variable:( ['time','lat','lon'], dat ) } )
		self.climatology.data = dat
		print( 'ds interpolated updated into self.ds' )
		return 1
	def downscale( self, output_dir, prefix=None ):
		import affine
		from affine import Affine
		import itertools
		from functools import partial
		from pathos.mp_map import mp_map

		operation_switch = { 'add':self.utils.add, 'mult':self.utils.mult }

		def two_digit_month( x ):
			''' make 1 digit month a standard 2-digit for output filenames '''
			month = str( x )
			if len(month) == 1:
				month = '0'+month
			return month

		time_suffix = [ '_'.join([two_digit_month( t.month ), str(t.year)]) for t in self.anomalies.time.to_pandas() ]

		# handle missing variable / model names
		if self.varname != None:
			variable = self.varname
		elif self.historical.variable != None:
			variable = self.historical.variable
		else:
			variable = 'variable'

		if self.modelname != None:
			model = self.modelname
		elif self.historical.model != None:
			model = self.historical.model
		else:
			model = 'model'

		output_filenames = [ os.path.join( output_dir, '_'.join([variable, self.historical.metric, self.historical.units, \
					self.historical.project, model, self.historical.scenario, ts]) + '.tif')  for ts in time_suffix ]

		# if there is a specific name prefix, use it
		if prefix != None:
			output_filenames = [ os.path.join( output_dir, '_'.join([prefix, ts]) + '.tif' ) for ts in time_suffix ]

		# rotate to pacific-centered
		if ( self.anomalies.lon.data > 200.0 ).any() == True:
			dat, lons = ( self.anomalies, self.anomalies.lon )
			self.anomalies_rot = dat
			src_transform = self.historical.transform_from_latlon( self.historical.ds.lat, lons )
			print( 'anomalies NOT rotated!' )
		else:
			dat, lons = self.utils.shiftgrid( 0., self.anomalies, self.anomalies.lon )
			self.anomalies_rot = dat
			src_transform = self.historical.transform_from_latlon( self.historical.ds.lat, lons )
			print( src_transform )
			print( 'anomalies rotated!' )

		# run and output
		rstlist = self.baseline.repeat( n=self.anomalies_rot.shape[0] / 12 ) # months
		
		if isinstance( self.anomalies_rot, xr.Dataset ):
			self.anomalies_rot = self.anomalies_rot[ self.historical.variable ].data
		elif isinstance( self.anomalies_rot, xr.DataArray ):
			self.anomalies_rot = self.anomalies_rot.data
		else:
			self.anomalies_rot = self.anomalies_rot

		args = zip( self.anomalies_rot, rstlist, output_filenames )

		args = [{'anom':i, 'base':j, 'output_filename':k,\
				'downscaling_operation':self.downscaling_operation, \
				'post_downscale_function':self.post_downscale_function,\
				'mask':self.mask, 'mask_value':self.mask_value } for i,j,k in args ]

		# partial and wrapper
		f = partial( self.utils.interp_ds, src_crs=self.src_crs, src_nodata=self.src_nodata, \
					dst_nodata=self.dst_nodata, src_transform=src_transform, resample_type=self.resample_type )

		run = partial( self.utils._run_ds, f=f, operation_switch=operation_switch, anom=self.anom, mask_value=self.mask_value )

		# run it
		out = mp_map( run, args, nproc=self.ncpus )
		return output_dir


# # # # # # # # # NEW FILL Dataset FOR A SPECIFIC SNAP ISSUE WITH pre DATA from CRU 
# # # # # # # # # # and pr DATA from CMIP5 / CRU
def find_boundary( arr ):
	'''
	return a mask of the boundary limit of DATA cells (overlays edge DATA not NA)
	this is especially useful if the data are land-only.  As in the CRU TS3.x data.
	'''
	from skimage.segmentation import find_boundaries
	bool_arr = np.copy( arr )
	ind = np.isnan( bool_arr )
	bool_arr[ ~ind ] = 1
	bool_arr[ ind ] = 0
	return find_boundaries( bool_arr, mode='inner' )

def calc_percentile( arr, aoi_mask=None, percentile=95, fill_value=0, nodata=None ):
	''' 
	calculate the percentile value potentially over a masked domain, 
	and avoiding nodata and np.nan AND return the nearest actual value to
	the np.nanpercentile( arr, percentile )

	arr = [numpy.ndarray] 2D array
	aoi_mask = [numpy.ndarray] 2D mask array of 0 (nomask) or 1 (mask)

	'''
	if aoi_mask is not None:
		# mask the background
		arr = arr[ (aoi_mask == fill_value) ]

	if nodata:
		arr = arr[ arr != nodata ]

	upperthresh = np.nanpercentile( arr, percentile )
	idx = (np.abs(arr - upperthresh)).argmin()
	return arr[ 1 ]

def correct_boundary( arr, bound_mask, aoi_mask=None, percentile=95, fill_value=0 ):
	''' correct the boundary pixels with non-acceptable values '''
	
	upperthresh = calc_percentile( arr, aoi_mask, 95, 0 )

	# drop any masks
	arr = np.array( arr )

	ind = np.where( bound_mask == True )
	vals = arr[ ind ]
	vals[ vals < 0.5 ] = 0.5
	vals[ vals > upperthresh ] = upperthresh
	arr[ ind ] = vals
	return np.array( arr )

def correct_inner( arr, bound_mask, aoi_mask=None, percentile=95, fill_value=0 ):
	''' correct the inner pixels with non-acceptable values '''

	upperthresh = calc_percentile( arr, aoi_mask, 95, 0 )
	
	# drop any masks
	arr = np.array( arr )

	ind = np.where( (arr > 0) & bound_mask != True )
	vals = arr[ ind ]
	vals[ vals < 0.5 ] = np.nan # set to the out-of-bounds value
	vals[ vals > upperthresh ] = upperthresh
	arr[ ind ] = vals
	return np.array( arr ) 

def correct_values( arr, aoi_mask=None, percentile=95, fill_value=0 ):
	''' correct the values for precip -- from @leonawicz'''

	upperthresh = calc_percentile( arr, aoi_mask, 95, 0 )

	# drop any masks
	arr = np.array( arr )

	arr[ arr < 0.5 ] = np.nan # set to the out-of-bounds value
	arr[ arr > upperthresh ] = upperthresh
	return np.array( arr ) 
