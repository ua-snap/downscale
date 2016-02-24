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
		metric='mean', ds_type='absolute', level=None, level_name=None ):
		'''
		this class appears to be working correctly for the CRU data as well
		this is great and means it should be built as a generic downscale class
		# a class for the AR5 would need to have
		1. historical
		2. future
		3. BaselineArr
		4. metric? --
		5. delta downscaletype
		6. 4th dimension -- plev?
		7. climatology begin / end
		'''
		self.historical = historical
		self.future = future
		self.baseline = baseline
		self.clim_begin = clim_begin
		self.clim_end = clim_end
		self.metric = metric
		self.ds_type = ds_type
		self.level = level
		self.level_name = level_name
		self.affine = self._calc_ar5_affine()
		self._concat_nc()
		self._calc_climatolgy()
		self._calc_anomalies()

	def _calc_ar5_affine( self, *args, **kwargs ):
		'''
		this assumes 0-360 longitude-ordering (pacific-centered)
		and WGS84 LatLong (Decimal Degrees). EPSG:4326.
		'''
		import affine
		lat_shape, lon_shape = self.historical.ds.dims[ 'lat' ], self.historical.ds.dims[ 'lon' ]
		lat_res = 180.0 / lat_shape
		lon_res = 360.0 / lon_shape
		return affine.Affine( lon_res, 0.0, 0.0, 0.0, -lat_res, 360.0 )
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
		if self.ds_type == 'absolute':
			anomalies = self.ds.groupby( 'time.month' ) - self.climatology
		elif self.ds_type == 'relative':
			anomalies = self.ds.groupby( 'time.month' ) / self.climatology
		else:
			NameError( '_calc_anomalies (ar5): value of ds_type must be "absolute" or "relative" ' )
		self.anomalies = anomalies
	@staticmethod
	def interp_ds( anom, base, output_filename, src_transform, downscaling_operation, post_downscale_function=None, mask=None, mask_value=0 ):
		'''	
		anom = [numpy.ndarray] 2-d array representing a single monthly timestep of the data to be downscaled. Must also be representative of anomalies.
		base = [str] filename of the corresponding baseline monthly file to use as template and downscale baseline for combining with anomalies.
		output_filename = [str] path to the output file to be created following downscaling
		src_transform = [affine.affine] 6 element affine transform of the input anomalies. [should be greenwich-centered]
		downscaling_operation = [str] one of 'add' or 'mult' depending on absolute or relative delta downscaling.
		post_downscale_function = [function] function that takes as input a single 2-d array and returns a 2-d array in the same shape as the input.  
		mask = [numpy.ndarray] 2-d array showing what values should be masked. Masked=0, unmasked=1. must be same shape as base.
		mask_value = [int] what value to use as the masked values. default=0.

		'''		
		from rasterio.warp import reproject, RESAMPLING
		# reproject / resample
		src_crs = {'init':'epsg:4326'}
		src_nodata = None # DangerTown™
		base = rasterio.open( base )
		baseline_arr = base.read( 1 )
		baseline_meta = base.meta
		baseline_meta.update( compress='lzw' )
		output_arr = np.empty_like( baseline_arr )

		# TODO: make this function available for manipulation if used for different needs
		reproject( anom, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
				dst_transform=baseline_meta['affine'], dst_crs=baseline_meta['crs'],\
				dst_nodata=None, resampling=RESAMPLING.cubic_spline, SOURCE_EXTRA=1000 )
		# downscale
		return utils.downscale( output_arr, baseline_arr, output_filename, downscaling_operation, \
				baseline_meta, post_downscale_function, mask, mask_value )
	def downscale( self, output_dir, prefix=None, ncpus=32 ):
		import affine
		import itertools
		from functools import partial
		# import multiprocessing
		from pathos import multiprocessing
		
		# output_filenames 
		time = self.anomalies.time.to_pandas()
		time_suffix = [ '_'.join([str(t.month), str(t.year)]) for t in time ]
		if prefix:
			output_filenames = [ os.path.join( output_dir, '_'.join([prefix, ts]) + '.tif' ) for ts in time_suffix ]
		else:
			if not self.historical.units:
				units = 'units'
			else:
				units = self.historical.units
			output_filenames = [ os.path.join( output_dir, '_'.join([self.historical.variable, units, \
							self.metric, self.historical.model, ts]) + '.tif')  for ts in time_suffix ]
		# rotate
		if ( self.anomalies.lon > 200.0 ).any() == True:
			dat, lons = utils.shiftgrid( 180., self.anomalies, self.anomalies.lon, start=False )
			a,b,c,d,e,f,g,h,i = self.affine 
			# flip it to the greenwich-centering
			src_transform = affine.Affine( a, b, -180.0, d, e, 180.0 )
		else:
			dat, lons = ( self.anomalies, self.anomalies.lon )
			src_transform = self.affine
		# run and output
		rstlist = self.baseline.filelist * (self.anomalies.shape[0] / 12)
		args = zip( self.anomalies, rstlist, output_filenames )
		args = [{'anom':i.data, 'base':j, 'output_filename':k} for i,j,k in args ]

		# determine operation type
		downscaling_operation_switch = {'absolute':'add', 'relative':'mult'}
		downscaling_operation = downscaling_operation_switch[ self.ds_type ]
		# partial and wrapper
		f = partial( self.interp_ds, src_transform=src_transform, downscaling_operation=downscaling_operation, \
						post_downscale_function=None, mask=None, mask_value=0 )
		def wrap( d ):
			return f( **d )

		pool = multiprocessing.Pool( ncpus )
		# out = pool.map( lambda x: f( **x ), args )
		out = pool.map( wrap, args )
		# pool.join()
		pool.close()
		return output_dir