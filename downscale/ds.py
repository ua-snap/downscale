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
	def __init__( self, baseline, clim_begin, clim_end, historical, future=None, \
		downscaling_operation='add', level=None, level_name=None, mask=None, mask_value=0, \
		ncpus=32, src_crs={'init':'epsg:4326'}, src_nodata=-9999.0, dst_nodata=None,
		post_downscale_function=None, varname=None, modelname=None, *args, **kwargs ):
		
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
		self.src_crs = src_crs
		self.src_nodata = src_nodata
		self.dst_nodata = dst_nodata
		self.post_downscale_function = post_downscale_function
		
		# calculate args
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
		if self.historical and self.future:
			self.anomalies = anomalies.sel( time=self.future.ds.time )
		else:
			self.anomalies = anomalies.sel( time=self.historical.ds.time )
	@staticmethod
	def interp_ds( anom, base, src_crs, src_nodata, dst_nodata, src_transform, *args, **kwargs ):
		'''	
		anom = [numpy.ndarray] 2-d array representing a single monthly timestep of the data to be downscaled. 
								Must also be representative of anomalies.
		base = [str] filename of the corresponding baseline monthly file to use as template and downscale 
								baseline for combining with anomalies.
		src_transform = [affine.affine] 6 element affine transform of the input anomalies. [should be greenwich-centered]

		'''		
		from rasterio.warp import reproject, RESAMPLING
		
		base = rasterio.open( base )
		baseline_arr = base.read( 1 )
		baseline_meta = base.meta
		baseline_meta.update( compress='lzw' )
		output_arr = np.empty_like( baseline_arr )

		reproject( anom, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
				dst_transform=baseline_meta['affine'], dst_crs=baseline_meta['crs'],\
				dst_nodata=dst_nodata, resampling=RESAMPLING.bilinear, SOURCE_EXTRA=5000 )
		return output_arr
	@staticmethod
	def wrap( d, f, operation_switch, anom=False, mask_value=0 ):
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

		# # # ANOMALIES TO DISK
		if anom == True:
			# write out the anomalies for testing:
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
		# # # ANOMALIES END

		dirname = os.path.dirname( d['output_filename'] )
		if not os.path.exists( dirname ):
			os.makedirs( dirname )

		# nested dict switch here, not attractive, but effective.
		output_arr = operation_switch[ d[ 'downscaling_operation' ] ]( base_arr, interped )
		
		# post downscale it if necessary:
		if post_downscale_function != None:
			output_arr = post_downscale_function( output_arr )
			# drop the mask if there is one -- this could become smarter
			if hasattr( output_arr, 'mask'):
				output_arr = output_arr.data

		# make sure data is masked
		output_arr[ mask == mask_value ] = meta[ 'nodata' ]

		# write it to disk.
		with rasterio.open( d[ 'output_filename' ], 'w', **meta ) as out:
			out.write( output_arr, 1 )
		return d['output_filename']
	def downscale( self, output_dir, prefix=None ):
		import affine
		import itertools
		from functools import partial
		from pathos.mp_map import mp_map
		
		# determine operation type
		def add( base, anom ):
			return base + anom
		def mult( base, anom ):
			return base * anom
		def div( base, anom ):
			# this one may not be useful, but the placeholder is here
			# return base / anom
			return NotImplementedError

		operation_switch = { 'add':add, 'mult':mult, 'div':div }

		def two_digit_month( x ):
			''' make 1 digit month a standard 2-digit for output filenames '''
			month = str( x )
			if len(month) == 1:
				month = '0'+month
			return month

		time_suffix = [ '_'.join([two_digit_month( t.month ), str(t.year)]) for t in self.anomalies.time.to_pandas() ]
		
		# # # # SOME STUFF TO DEAL WITH CUSTOM OUTPUT NAMING OTF # # # # 
		# deal with missing variable names and/or model names
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

		# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
		# set up some output filenames
		output_filenames = [ os.path.join( output_dir, '_'.join([variable, self.historical.metric, self.historical.units, \
					self.historical.project, model, self.historical.scenario, ts]) + '.tif')  for ts in time_suffix ]

		# if there is a specific name prefix, use it
		if prefix != None:
			output_filenames = [ os.path.join( output_dir, '_'.join([prefix, ts]) + '.tif' ) for ts in time_suffix ]
		
		# rotate to greenwich-centered
		if ( self.anomalies.lon.data > 200.0 ).any() == True:
			dat, lons = utils.shiftgrid( 180., self.anomalies, self.anomalies.lon, start=False )
			self.anomalies_rot = dat
			a,b,c,d,e,f,g,h,i = self.affine
			# flip it to the greenwich-centering
			src_transform = affine.Affine( a, b, -180.0, d, e, 90.0 )
			print( 'anomalies rotated!' )
		else:
			dat, lons = ( self.anomalies, self.anomalies.lon )
			self.anomalies_rot = dat 
			src_transform = self.affine
			print( 'anomalies NOT rotated!' )
	
		# run and output
		rstlist = self.baseline.filelist * (self.anomalies_rot.shape[0] / 12)
		args = zip( self.anomalies_rot, rstlist, output_filenames )

		args = [{'anom':i, 'base':j, 'output_filename':k,\
				'downscaling_operation':self.downscaling_operation, \
				'post_downscale_function':self.post_downscale_function,\
				'mask':self.mask, 'mask_value':self.mask_value } for i,j,k in args ]

		# partial and wrapper
		f = partial( self.interp_ds, src_crs=self.src_crs, src_nodata=self.src_nodata, \
					dst_nodata=self.dst_nodata, src_transform=src_transform )
		
		wrapped = partial( self.wrap, f=f, operation_switch=operation_switch, anom=True, mask_value=self.mask_value )
		
		# run it
		out = mp_map( wrapped, args, nproc=self.ncpus )
		return output_dir
