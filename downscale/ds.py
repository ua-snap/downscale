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
import xarray as xr
from downscale import utils

class DeltaDownscale( object ):
	def __init__( self, baseline, clim_begin, clim_end, historical, future=None,
				downscaling_operation='add', level=None, level_name=None, 
				mask=None, mask_value=0,ncpus=32, src_crs={'init':'epsg:4326'}, 
				src_nodata=-9999.0, dst_nodata=None, post_downscale_function=None, 
				varname=None, modelname=None, anom=False, resample_type='bilinear', 
				fix_clim=False, interp=False, *args, **kwargs ):
		
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

		# calculate args
		self.anomalies = None
		self.ds = None
		self._concat_nc()
		self._calc_climatolgy()
		# fix pr climatologies if desired
		if fix_clim == True: # ARG ME!
			self._fix_clim()

		# calculate anomalies with the new climatology values
		self._calc_anomalies()

		# interpolate across space here instead of in `Dataset`
		self._rotated = False # brought from dataset KEEP?
		self._lonpc = None # brought from dataset KEEP?
		if interp == True or fix_clim == True:
			print( 'running interpolation across NAs -- base resolution' )
			_ = self.interp_na( )

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

	# # # # # # # # # MOVED FROM Dataset
	def _fix_clim( self, find_bounds=False ):
		''' fix values in precip data '''
		if find_bounds == True:
			bound_mask = find_boundary( climatology[ 0, ... ].data )
			for idx in range( self.climatology.shape[0] ):
				arr = self.climatology[ idx, ... ].data
				arr = correct_boundary( arr, bound_mask )
				self.climatology[ idx, ... ].data = correct_inner( arr, bound_mask )

		elif find_bounds == False:
			for idx in range( self.climatology.shape[0] ):
				arr = self.climatology[ idx, ... ].data
				self.climatology[ idx, ... ].data = correct_values( arr )
		else:
			ValueError( 'find_bounds arg is boolean only' )
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
			dat, lons = self.ds[ self.historical.variable ].data, self.ds.lon
			self._lonpc = lons
		else:
			# greenwich-centered rotate to 0-360 for interpolation across pacific
			dat, lons = self.rotate( self.ds[ self.historical.variable ].values, self.ds.lon, to_pacific=True )
			self._rotated = True # update the rotated attribute
			self._lonpc = lons

		# mesh the lons and lats and unravel them to 1-D
		xi, yi = np.meshgrid( self._lonpc, self.ds.lat.data )
		lo, la = [ i.ravel() for i in (xi,yi) ]

		# setup args for multiprocessing
		df_list = [ pd.DataFrame({ 'x':lo, 'y':la, 'z':d.ravel() }).dropna( axis=0, how='any' ) for d in dat ]

		args = [ {'x':np.array(df['x']), 'y':np.array(df['y']), 'z':np.array(df['z']), \
				'grid':(xi,yi), 'method':self.historical.method, 'output_dtype':output_dtype } for df in df_list ]
		
		# # # # USE MLAB's griddata which we _can_ parallelize
		def wrap( d ):
			''' simple wrapper around utils.xyz_to_grid for mp_map '''
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
		self.ds = self.ds.update( { self.historical.variable:( ['time','lat','lon'], dat ) } )
		print( 'ds interpolated updated into self.ds' )
		return 1
	# # # # # # # # # END! MOVED FROM Dataset
	@staticmethod
	def interp_ds( anom, base, src_crs, src_nodata, dst_nodata, src_transform, resample_type='bilinear',*args, **kwargs ):
		'''	
		anom = [numpy.ndarray] 2-d array representing a single monthly timestep of the data to be downscaled. 
								Must also be representative of anomalies.
		base = [str] filename of the corresponding baseline monthly file to use as template and downscale 
								baseline for combining with anomalies.
		src_transform = [affine.affine] 6 element affine transform of the input anomalies. [should be greenwich-centered]
		resample_type = [str] one of ['bilinear', 'count', 'nearest', 'mode', 'cubic', 'index', 'average', 'lanczos', 'cubic_spline']
		'''	
		import rasterio
		from rasterio.warp import reproject, RESAMPLING

		resampling = {'average':RESAMPLING.average,
					'cubic':RESAMPLING.cubic,
					'lanczos':RESAMPLING.lanczos,
					'bilinear':RESAMPLING.bilinear,
					'cubic_spline':RESAMPLING.cubic_spline,
					'mode':RESAMPLING.mode,
					'count':RESAMPLING.count,
					'index':RESAMPLING.index,
					'nearest':RESAMPLING.nearest }
		
		base = rasterio.open( base )
		baseline_arr = base.read( 1 )
		baseline_meta = base.meta
		baseline_meta.update( compress='lzw' )
		output_arr = np.empty_like( baseline_arr )
		
		reproject( anom, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
				dst_transform=baseline_meta['affine'], dst_crs=baseline_meta['crs'],\
				dst_nodata=dst_nodata, resampling=resampling[ resample_type ], SOURCE_EXTRA=1000 )
		return output_arr
	@staticmethod
	def _run_ds( d, f, operation_switch, anom=False, mask_value=0 ):
		'''
		run the meat of downscaling with this runner function for parallel processing

		ARGUMENTS:
		----------
		d = [dict] kwargs dict of args to pass to interpolation function
		f = [ ]

		RETURNS:
		--------

		'''
		import copy
		post_downscale_function = d[ 'post_downscale_function' ]
		interped = f( **d )
		base = rasterio.open( d[ 'base' ] )
		base_arr = base.read( 1 )
		mask = base.read_masks( 1 )

		# set up output file metadata.
		meta = base.meta
		meta.update( compress='lzw' )
		if 'transform' in meta.keys():
			meta.pop( 'transform' )

		if anom == True:
			# write out the anomalies
			anom_filename = copy.copy( d[ 'output_filename' ] )
			dirname, basename = os.path.split( anom_filename )
			dirname = os.path.join( dirname, 'anom' )
			basename = basename.replace( '.tif', '_anom.tif' )
			try:
				if not os.path.exists( dirname ):
					os.makedirs( dirname )
			except:
				pass
			anom_filename = os.path.join( dirname, basename )
			with rasterio.open( anom_filename, 'w', **meta ) as anom:
				anom.write( interped, 1 )
		
		# make sure the output dir exists and if not, create it
		dirname = os.path.dirname( d[ 'output_filename' ] )
		if not os.path.exists( dirname ):
			os.makedirs( dirname )

		# operation switch
		output_arr = operation_switch[ d[ 'downscaling_operation' ] ]( base_arr, interped )
		
		# post downscale it if func given
		if post_downscale_function != None:
			output_arr = post_downscale_function( output_arr )
			# drop the mask if there is one
			if hasattr( output_arr, 'mask'):
				output_arr = output_arr.data

		# make sure data is masked
		output_arr[ mask == mask_value ] = meta[ 'nodata' ]

		# write it to disk.
		with rasterio.open( d[ 'output_filename' ], 'w', **meta ) as out:
			out.write( output_arr, 1 )
		return d['output_filename']
	@staticmethod
	def add( base, anom ):
		''' add anomalies to baseline '''
		return base + anom
	@staticmethod
	def mult( base, anom ):
		''' multiply anomalies to baseline '''
		return base * anom
	def downscale( self, output_dir, prefix=None ):
		import affine
		from affine import Affine
		import itertools
		from functools import partial
		from pathos.mp_map import mp_map

		operation_switch = { 'add':self.add, 'mult':self.mult }

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
			dat, lons = utils.shiftgrid( 0., self.anomalies, self.anomalies.lon )
			self.anomalies_rot = dat
			src_transform = self.historical.transform_from_latlon( self.historical.ds.lat, lons )
			print( src_transform )
			print( 'anomalies rotated!' )

		# run and output
		rstlist = self.baseline.filelist * (self.anomalies_rot.shape[0] / 12)
		
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
		f = partial( self.interp_ds, src_crs=self.src_crs, src_nodata=self.src_nodata, \
					dst_nodata=self.dst_nodata, src_transform=src_transform, resample_type=self.resample_type )

		run = partial( self._run_ds, f=f, operation_switch=operation_switch, anom=self.anom, mask_value=self.mask_value )

		# run it
		out = mp_map( run, args, nproc=self.ncpus )
		return output_dir


# # # # # # # # # NEW FILL Dataset FOR A SPECIFIC SNAP ISSUE WITH pre DATA from CRU 
# # # # # # # # # # and pr DATA from CMIP5
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

def correct_boundary( arr, bound_mask, percentile=95 ):
	''' correct the boundary pixels with non-acceptable values '''
	upperthresh = np.percentile( arr[~np.isnan( arr )], percentile )
	ind = np.where( bound_mask == True )
	vals = arr[ ind ]
	vals[ vals < 0.5 ] = 0.5
	vals[ vals > upperthresh ] = upperthresh
	arr[ ind ] = vals
	return arr

def correct_inner( arr, bound_mask, percentile=95 ):
	''' correct the inner pixels with non-acceptable values '''
	upperthresh = np.percentile( arr[~np.isnan( arr )], percentile )
	mask = np.copy( arr )	
	ind = np.where( (arr > 0) & bound_mask != True )
	vals = arr[ ind ]
	vals[ vals < 0.5 ] = np.nan # set to the out-of-bounds value
	vals[ vals > upperthresh ] = upperthresh
	arr[ ind ] = vals
	return arr

def correct_values( arr, percentile=95 ):
	''' correct the values for precip -- from @leonawicz'''
	upperthresh = np.percentile( arr[~np.isnan( arr )], percentile )
	arr[ arr < 0.5 ] = np.nan # set to the out-of-bounds value
	arr[ arr > upperthresh ] = upperthresh
	return arr

# # # # # # # # # END! NEW FILL Dataset

