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
import xarray as xr
from downscale import utils

class DeltaDownscale( object ):
	def __init__( self, baseline, clim_begin, clim_end, historical, future=None, \
		downscaling_operation='add', level=None, level_name=None, mask=None, mask_value=0, \
		ncpus=32, src_crs={'init':'epsg:4326'}, src_nodata=-9999.0, dst_nodata=None,
		post_downscale_function=None, varname=None, modelname=None ):
		
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
		self.affine = self.historical._calc_affine()

		# new args
		self.src_crs = src_crs
		self.src_nodata = src_nodata
		self.dst_nodata = dst_nodata
		self.post_downscale_function = post_downscale_function
		
		self.anomalies = None
		self._concat_nc()
		self._calc_climatolgy()
		self._calc_anomalies()

	def _concat_nc( self ):
		if self.historical and self.future:
			ds = xr.concat([ self.historical.ds, self.future.ds ], dim='time' )
		else:
			ds = self.historical.ds
		ds = ds[ self.historical.variable ]
		if self.level:
			levidx, = np.where( ds[ self.level_name ] == self.level )
			ds = ds[ :, levidx[0], ... ]
		self.ds = ds
	def _calc_climatolgy( self ):
		'''slice / aggregate to climatology using mean'''
		try:
			climatology = self.ds.loc[ { 'time' : slice( self.clim_begin, self.clim_end ) } ]
			self.climatology = climatology.groupby( 'time.month' ).mean( 'time' )
		except Exception:
			raise AttributeError( 'non-overlapping climatology period and series' )
	def _calc_anomalies( self ):
		# anomalies
		if self.downscaling_operation == 'add':
			anomalies = self.ds.groupby( 'time.month' ) - self.climatology
		elif self.downscaling_operation == 'mult':
			anomalies = self.ds.groupby( 'time.month' ) / self.climatology
		else:
			NameError( '_calc_anomalies (ar5): value of downscaling_operation must be "add" or "mult" ' )
		self.anomalies = anomalies
	@staticmethod
	def interp_ds( anom, base, src_crs, src_nodata, dst_nodata, src_transform, *args, **kwargs ):
		'''	
		anom = [numpy.ndarray] 2-d array representing a single monthly timestep of the data to be downscaled. Must also be representative of anomalies.
		base = [str] filename of the corresponding baseline monthly file to use as template and downscale baseline for combining with anomalies.
		# REMOVE output_filename = [str] path to the output file to be created following downscaling
		src_transform = [affine.affine] 6 element affine transform of the input anomalies. [should be greenwich-centered]
		# REMOVE downscaling_operation = [str] one of 'add' or 'mult' depending on absolute or relative delta downscaling.
		# REMOVE post_downscale_function = [function] function that takes as input a single 2-d array and returns a 2-d array in the same shape as the input.  
		# REMOVE mask = [numpy.ndarray] 2-d array showing what values should be masked. Masked=0, unmasked=1. must be same shape as base.
		# REMOVE mask_value = [int] what value to use as the masked values. default=0.

		'''		
		from rasterio.warp import reproject, RESAMPLING
		# reproject / resample
		base = rasterio.open( base )
		baseline_arr = base.read( 1 )
		baseline_meta = base.meta
		baseline_meta.update( compress='lzw' )
		output_arr = np.empty_like( baseline_arr )

		reproject( anom, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
				dst_transform=baseline_meta['affine'], dst_crs=baseline_meta['crs'],\
				dst_nodata=dst_nodata, resampling=RESAMPLING.cubic_spline, SOURCE_EXTRA=1000 )
		return output_arr
	def downscale( self, output_dir, prefix=None ):
		import affine
		import itertools
		from functools import partial
		from pathos import multiprocessing
		
				# determine operation type
		def add( base, anom ):
			return base + anom
		def mult( base, anom ):
			return base * anom
		def div( base, anom ):
			# this one may not be useful, but the placeholder is here
			# return base / anom
			return NotImplementedError
		try:
			operation_switch = { 'add':add, 'mult':mult, 'div':div }
		except:
			AttributeError( 'downscale: incorrect downscaling_operation str' )

		time_arr = self.anomalies.time.to_pandas() # CHANGE THIS NAME!

		# we need to be able to output ONLY the years we want if there is a future
		if self.future:
			time_arr = self.future.ds.time.to_pandas()
		
		# slice the anomalies
		self.anomalies = self.anomalies.sel( time=slice( time_arr[0], time_arr[-1] ) )
		
		time_suffix = [ '_'.join([str(t.month), str(t.year)]) for t in time_arr ]
		
		# deal with missing variable names and/or model names
		if self.historical.variable != None:
			variable = self.historical.variable
		else:
			variable = 'variable'
		
		if self.historical.model != None:
			model = self.historical.model
		else:
			model = 'model'

		if self.modelname != None:
			model = self.modelname

		# set up some output filenames
		# # # [ variable, metric, units, project, model, scenario, month, year, ext ]
		output_filenames = [ os.path.join( output_dir, '_'.join([variable, self.historical.metric, self.historical.units, \
					self.historical.project, model, self.historical.scenario, ts]) + '.tif')  for ts in time_suffix ]

		# if there is a specific name prefix, use it
		if prefix != None:
			output_filenames = [ os.path.join( output_dir, '_'.join([prefix, ts]) + '.tif' ) for ts in time_suffix ]
		
		# rotate
		if ( self.anomalies.lon.data > 200.0 ).any() == True:
			dat, lons = utils.shiftgrid( 180., self.anomalies, self.anomalies.lon, start=False )
			a,b,c,d,e,f,g,h,i = self.affine
			# flip it to the greenwich-centering
			src_transform = affine.Affine( a, b, -180.0, d, e, 90.0 )
			print( 'anomalies rotated!' )
		else:
			dat, lons = ( self.anomalies, self.anomalies.lon )
			src_transform = self.affine
			print( 'anomalies NOT rotated!' )
	
		# run and output
		rstlist = self.baseline.filelist * (self.anomalies.shape[0] / 12)
		args = zip( self.anomalies, rstlist, output_filenames )

		args = [{'anom':i.data, 'base':j, 'output_filename':k,\
				'downscaling_operation':self.downscaling_operation, \
				'post_downscale_function':self.post_downscale_function,\
				'mask':self.mask, 'mask_value':self.mask_value } for i,j,k in args ]

		# downscaling_operation = operation_switch[ self.downscaling_operation ]
		# partial and wrapper
		f = partial( self.interp_ds, src_crs=self.src_crs, src_nodata=self.src_nodata, dst_nodata=self.dst_nodata, \
					src_transform=src_transform )

		def wrap( d ):
			interped = f( **d )
			base = rasterio.open( d[ 'base' ] )
			base_arr = base.read( 1 )
			tmp = base.read_masks( 1 )
			output_arr = operation_switch[ d[ 'downscaling_operation' ] ]( base_arr, interped )

			# mask
			output_arr[ tmp == 0 ] = base.nodata

			# output_arr[ np.isinf( output_arr ) ] = meta[ 'nodata' ] # not sure about this one

			# post downscale it if necessary:
			if d['post_downscale_function'] != None:
				output_arr = post_downscale_function( output_arr )

			meta = base.meta
			meta.update( compress='lzw' )
			if 'transform' in meta.keys():
				meta.pop( 'transform' )

			with rasterio.open( d[ 'output_filename' ], 'w', **meta ) as out:
				out.write( output_arr, 1 )
			return d['output_filename']

		pool = multiprocessing.Pool( self.ncpus )
		out = pool.map( wrap, args )
		pool.close()
		pool.join()
		return output_dir
