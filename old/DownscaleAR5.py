# -*- coding: utf8 -*-
# # #
# Downscale PCMDI AR5 data to a pre-processed climatology
#  extent, resolution, reference system
#
# Author: Michael Lindgren (malindgren@alaska.edu)
# # #
import rasterio, os
import numpy as np
from downscale import DownscalingUtils as utils #utils

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
		self.filelist = filelist
		self.meta = rasterio.open( self.filelist[0] ).meta
		self.arrlist = [ rasterio.open( fn ).read( 1 ) for fn in self.filelist ]

class Dataset( object ):
	def __init__( self, fn, variable, model, scenario, units=None, *args, **kwargs ):
		'''
		fn = [str] path to the xray supported dataset to be read in.
		variable = [str] abbreviation of variable name to extract from file
		model = [str] name of the model being read
		scenario = [str] name of the scenario being read
		units = [str] abbreviation of the units of the variable
		'''
		import xarray as xr
		self.fn = fn
		self.ds = xr.open_dataset( fn )
		self.variable = variable
		self.model = model
		self.scenario = scenario
		self.units = units

class AR5( object ):
	def __init__( self, baseline, clim_begin, clim_end, historical, future=None, \
		metric='mean', ds_type='absolute', level=None, level_name=None ):
		'''	
		this class appears to be working correctly for the CRU data as well
		this is great and means it should be built as a generic downscale class
		# a class for the AR5 would need to have
		1. historical
		2. modeled
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
		import xarray as xr
		import affine
		lat_shape, lon_shape = self.historical.ds.dims[ 'lat' ], self.historical.ds.dims[ 'lon' ]
		lat_res = 180.0 / lat_shape
		lon_res = 360.0 / lon_shape
		return affine.Affine( lon_res, 0.0, 0.0, 0.0, -lat_res, 360.0 )
	def _concat_nc( self ):
		if self.historical and self.future:
			ds = xr.concat([ self.historical.ds, self.modeled.ds ])
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
		from rasterio.warp import reproject, RESAMPLING
		# reproject / resample
		src_crs = {'init':'epsg:4326'}
		src_nodata = None # DangerTown™
		baseline_meta = base.meta
		baseline_meta.update( compress='lzw' )
		baseline_arr = base.read( 1 )
		output_arr = np.empty_like( base.read( 1 ) ) # base should be a rasterio obj not arr... check this.

		# TODO: make this function available for manipulation if used for different needs
		reproject( output_arr, baseline_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
				dst_transform=baseline_meta['affine'], dst_crs=baseline_meta['crs'],\
				dst_nodata=None, resampling=RESAMPLING.cubic_spline, SOURCE_EXTRA=1000 )
		# downscale
		return utils.downscale( anom, output_arr, output_filename, downscaling_operation, \
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
		rstlist = [ rasterio.open( fn ) for fn in self.baseline.filelist ]
		# # # # # 
		# maybe itertools.repeat is needed here? # these are same references
		rstlist = itertools.repeat( rstlist, (self.anomalies.shape[0] / 12) )
		rstlist = [ j for i in rstlist for j in i ]
		# rstlist = rstlist * (self.anomalies.shape[0] / 12)
		# # # # # 
		args = zip( self.anomalies, rstlist, output_filenames )
		args = [{'anom':i.data, 'base':j, 'output_filename':k} for i,j,k in args ]

		downscaling_operation_switch = {'absolute':'add', 'relative':'mult'}
		downscaling_operation = downscaling_operation_switch[ self.ds_type ]
		
		# anom_arr, baseline_arr, output_filename,	downscaling_operation, \
		# 	meta, post_downscale_function, mask=None, mask_value=0

		f = partial( self.interp_ds, src_transform=src_transform, downscaling_operation=downscaling_operation, post_downscale_function=None, mask=None, mask_value=0 )
		pool = multiprocessing.Pool( ncpus )
		out = pool.map( lambda x: f( **x ), args[:11] )
		pool.join()
		pool.close()
		return output_dir



#  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 




def downscale_wrapper( arr, affine, crs, baseline, output_filename, downscaling_operation, post_downscale_function ):
	# rotate
	if ( self.data.ds.lon > 200.0 ).any() == True:
		dat, lons = utils.shiftgrid( 180., self.data.anomalies, self.historical.ds.lon, start=False )
		a,b,c,d,e,f,g,h,i = affine #flip it to the greenwich-centering
		src_transform = affine.Affine( a, b, -180.0, d, e, 180.0 )
	else:
		dat, lons = ( self.data.ds, self.historical.ds.lon )
		src_transform = affine
	
	# reproject / resample
	src_crs = {'init':'epsg:4326'}
	src_nodata = None # DangerTown™
	baseline_meta = baseline.meta
	baseline_meta.update( compress='lzw' )
	output_arr = np.empty_like( baeline.read( 1 ) )

	# TODO: make this function available for manipulation if used for different needs
	reproject( arr, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
			dst_transform=baseline_meta['affine'], dst_crs=baseline_meta['crs'],\
			dst_nodata=None, resampling=RESAMPLING.cubic_spline, SOURCE_EXTRA=1000 )
	
	# downscale
	return utils.downscale( arr, output_arr, output_filename, downscaling_operation, \
					baseline_meta, post_downscale_function, mask=None, mask_value=0 )


def run_downscale( data, baseline, output_path, postprep_function, ncpus=2 ):
	from functools import partial
	output_filenames = [ os.path.join( output_path, monyear + '.tif' ) for monyear in monyear ]
	data_base_fn = zip( data.ds, baseline.arr_list, output_filenames )
	pool = multiprocessing.Pool( ncpus )
	f = partial( downscale_wrapper,  )
	out = pool.map( f, data_base_fn )




# def downscale( data, baseline, output_filenames, postprep_function, ncpus=2 ):
# 	from downscale import utils
# 	import multiprocessing





# class Downscale( object ):
# 	''' downscaling... '''
# 	def __init__( self, data, baseline, *args, **kwargs ):
# 		self.data = data
# 		self.baseline = baseline
# 		self.historical = self.data.historical
# 	def rotate( self ):
# 		'''rotate globe back to -180.0 to 180.0 longitudes'''
# 		if ( self.historical.ds.lon > 200.0 ).any() == True:
# 			dat, lons = utils.shiftgrid( 180., self.data.anomalies, self.historical.ds.lon, start=False )
# 		else:
# 			dat, lons = ( self.historical.ds.data, self.historical.ds.lon )
# 		return dat, lons
# 	def reproject( self ):
# 		'''	reproject / resample / crop to baseline	'''
# 		rst = rasterio.open( self.baseline.filelist[0] )
# 		output_arr = np.empty_like( rst.read( 1 ) )
# 		a,b,c,d,e,f,g,h,i = self.affine #flip it to the greenwich-centering
# 		src_transform = affine.Affine( a, b, -180.0, d, e, 180.0 )
# 		src_crs = {'init':'epsg:3338'}
# 		src_nodata = None # DangerTown™
# 		baseline_meta = self.baseline.meta
# 		baseline_meta.update( compress='lzw' )
# 		# rotate if needed:
# 		dat, lons = self.rotate( )

# 		dat_list = [ arr for arr in dat ]

# 		# reproject it
# 		reproject( dat, output_arr, src_transform=self.affine, src_crs=src_crs, src_nodata=src_nodata, \
# 				dst_transform=baseline_meta['affine'], dst_crs=baseline_meta['crs'],\
# 				dst_nodata=None, resampling=RESAMPLING.cubic_spline, SOURCE_EXTRA=1000 )

# 		# mask it with the internal mask in the template raster, where 0 is oob. DangerTown™
# 		mask = rst.read_masks( 1 ) == 0
# 		output_arr[ mask ] = baseline_meta[ 'nodata' ]

# 	# now we need to use the downscale function from the utilities module and write to disk.




