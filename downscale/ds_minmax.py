# -*- coding: utf8 -*-
# # #
# Downscale PCMDI AR5 data to a pre-processed climatology
#  extent, resolution, reference system
#
#  ::NOTE:: this version of DeltaDownscale is built for tmin/tmax
#	data where our old method causes for the mins and max and means to 
# 	cross each other in non-normal ways.
#
# Author: Michael Lindgren (malindgren@alaska.edu)
# # #

from downscale import DeltaDownscale, utils
import os, rasterio
import numpy as np
import xarray as xr

def delta_mm( fn, mean_fn, variable, mean_variable='tas' ):
	'''
	simple way to compute extreme - mean deltas as 
	native model resolution and write to NetCDF4 on disk
	'''
	ds = xr.open_dataset( fn )[ variable ]
	ds_mean = xr.open_dataset( mean_fn )[ mean_variable ]
	delta = ds - ds_mean
	return delta.to_dataset( name=variable )

class DeltaDownscaleMinMax( DeltaDownscale ):
	def __init__( self, mean_ds=None, mean_variable=None, *args, **kwargs ): 
		'''
			note here that all data falls into the 'historical' category, because we no longer need to 
			have the 1961-1990 climatology period for the futures as this version of DeltaDownscale computes
			deltas by removing the mean in time instead of removing the climatology.
		'''
		# if there is no mean dataset to work with --> party's over
		print('mean:{}'.format(mean_ds))
		
		if mean_ds is None:
			raise Exception( 'you must include the mean variable in the raw resolution \
								as arg `mean_ds`=downscale.Dataset object or use `DeltaDownscale`' )
		# setup new args
		self.mean_ds = mean_ds.ds[ mean_variable ] # new
		self.mean_variable = mean_variable

		print( kwargs.keys() )

		# force a false for interpolation of NA's with Super...
		if 'interp' in kwargs.keys():
			print( 'kwargs-interpval:{}'.format(kwargs['interp']) )
			interp_val = kwargs.pop( 'interp' )
			kwargs.update( interp=False )
			print( 'interp_val:{}'.format(interp_val) )

		# remove the super call since this is python2 and it suuucks...
		# DeltaDownscale.__init__( self, **kwargs )

		super( DeltaDownscaleMinMax, self ).__init__( **kwargs )
		# print('finished super()!')

		if 'interp' in kwargs: 
			# now reset the interpolation value so we can interpolate the anomalies
			# INSTEAD of interpolating the input data series...  This may be better practice overall.
			# NOT changing it.
			self.interp = interp_val

		# mask some properties from the super() class. that are unneeded.
		self.anomalies = None
		self.clim_begin = None
		self.clim_end = None

		# calc deltas between the mean and the extreme data set 
		print( 'calc anoms minmax' )
		self._calc_anomalies()
		
		# # TESTING
		# print('type_mean_ds: {} '.format( type( self.mean_ds ) ) )
		print( 'self.interp: {}'.format(self.interp) )
		if self.interp == True:
			print( 'running interpolation across NAs -- base resolution -- !ANOMALIES! dataset' )
			self.interp_na( )

	def _calc_climatolgy( self ):
		''' MASK THIS FOR MINMAX slice / aggregate to climatology using mean'''
		self.climatology = None
	def _calc_anomalies( self ):
		''' calculate deltas but call them anomalies to fit the `downscale` pkg methods '''			
		if self.downscaling_operation == 'add':
			print( 'calc_anom minmax version' )
			# anomalies = (self.historical.ds[ self.historical.variable ] - self.mean_ds.ds[ self.mean_variable ] ) #.to_dataset( name=variable )
			self.anomalies = (self.ds - self.mean_ds ) #.to_dataset( name=variable )
		elif self.downscaling_operation == 'mult':
			# anomalies = (self.historical.ds[ self.historical.variable ] / self.mean_ds.ds[ self.mean_variable ] ) #.to_dataset( name=variable )
			self.anomalies = (self.ds / self.mean_ds ) #.to_dataset( name=variable )
		else:
			NameError( '_calc_anomalies (ar5): value of downscaling_operation must be "add" or "mult" ' )
		# self.mean_ds = None # watch this one... trying to save on RAM... 
	# def _interp_na_mean( self ):
	# 	'''
	# 	np.float32
	# 	method = [str] one of 'cubic', 'near', 'linear'

	# 	return a list of dicts to pass to the xyz_to_grid in parallel
	# 	'''
	# 	from copy import copy
	# 	import pandas as pd
	# 	import numpy as np
	# 	from pathos.mp_map import mp_map

	# 	# remove the darn scientific notation
	# 	np.set_printoptions( suppress=True )
	# 	output_dtype = np.float32
		
	# 	# if 0-360 leave it alone
	# 	if ( self.mean_ds.lon > 200.0 ).any() == True:
	# 		dat, lons = self.mean_ds.data, self.mean_ds.lon
	# 		self._lonpc = lons
	# 	else:
	# 		# greenwich-centered rotate to 0-360 for interpolation across pacific
	# 		dat, lons = self.utils.rotate( self.mean_ds.values, self.mean_ds.lon, to_pacific=True )
	# 		self._rotated = True # update the rotated attribute
	# 		self._lonpc = lons

	# 	# mesh the lons and lats and unravel them to 1-D
	# 	xi, yi = np.meshgrid( self._lonpc, self.mean_ds.lat.data )
	# 	lo, la = [ i.ravel() for i in (xi,yi) ]

	# 	# setup args for multiprocessing
	# 	df_list = [ pd.DataFrame({ 'x':lo, 'y':la, 'z':d.ravel() }).dropna( axis=0, how='any' ) for d in dat ]

	# 	args = [ {'x':np.array(df['x']), 'y':np.array(df['y']), 'z':np.array(df['z']), \
	# 			'grid':(xi,yi), 'method':self.historical.method, 'output_dtype':output_dtype } for df in df_list ]
		
	# 	print( 'processing interpolation to convex hull in parallel using {} cpus.'.format( self.ncpus ) )
	# 	dat_list = mp_map( self.wrap, args, nproc=self.ncpus )
	# 	dat_list = [ np.array(i) for i in dat_list ] # drop the output mask
	# 	dat = np.array( dat_list )

	# 	lons = self._lonpc
	# 	if self._rotated == True: # rotate it back
	# 		dat, lons = self.utils.rotate( dat, lons, to_pacific=False )
				
	# 	# place back into a new xarray.Dataset object for further processing
	# 	# self.mean_ds = self.mean_ds.update( { self.historical.variable:( ['time','lat','lon'], dat ) } )
	# 	self.mean_ds.data = dat
	# 	print( 'ds interpolated updated into self.mean_ds' )
	# 	return 1

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
		if ( np.array(self.anomalies.lon) > 200.0 ).any() == True:
			dat, lons = np.array(self.anomalies.data), np.array(self.anomalies.lon)
			self._lonpc = lons
		else:
			# greenwich-centered rotate to 0-360 for interpolation across pacific
			dat, lons = self.utils.rotate( np.array(self.anomalies.values), np.array(self.anomalies.lon), to_pacific=True )
			self._rotated = True # update the rotated attribute
			self._lonpc = lons

		# mesh the lons and lats and unravel them to 1-D
		xi, yi = np.meshgrid( self._lonpc, self.anomalies.lat.data )
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
		# self.anomalies = self.anomalies.update( { self.historical.variable:( ['time','lat','lon'], dat ) } )
		self.anomalies.data = dat
		print( 'anomalies interpolated updated into self.anomalies' )
		return 1	
	def downscale( self, output_dir, prefix=None ):
		'''
		updated version of downscale function to mask the non-minmax version and how
		it works with baseline climatology vs. the full mean series as with the min/max
		'''
		import affine, rasterio
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
			dat, lons = ( np.array(self.anomalies), np.array(self.anomalies.lon) )
			self.anomalies_rot = dat
			src_transform = self.historical.transform_from_latlon( self.ds.lat, lons )
			# print( 'anomalies NOT rotated!' )
		else:
			dat, lons = self.utils.shiftgrid( 0., np.array(self.anomalies), np.array(self.anomalies.lon) )
			self.anomalies_rot = dat
			src_transform = self.historical.transform_from_latlon( self.ds.lat, lons )
			# print( 'anomalies rotated!' )

		# # # # # #TSSTING STUFF
		# count, height, width = dat.shape
		# meta = { 'dtype':'float32', 'driver':'GTiff', 'count':1, 'width':width, 'height':height, 'compress':'lzw', 'affine':src_transform }
		# with rasterio.open( '/workspace/UA/malindgren/temporary/TEST_REGRID.tif', 'w', **meta ) as out:
		# 	out.write( dat[0].astype( np.float32 ), 1 )
		# # # # # # # # # # # # 

		# # # IMPORTANT: list all files since it without a REPEAT since it is tasmin/max...
		rstlist = self.baseline.filelist
		
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
		f = partial( self.utils.interp_ds, src_crs=self.src_crs, src_nodata=None, \
					dst_nodata=None, src_transform=src_transform, resample_type=self.resample_type )

		run = partial( self.utils._run_ds, f=f, operation_switch=operation_switch, anom=self.anom, mask_value=self.mask_value )

		# run it
		out = mp_map( run, args, nproc=self.ncpus )
		return output_dir
	# @staticmethod
	# def interp_ds( anom, base, src_crs, src_nodata, dst_nodata, src_transform, resample_type='bilinear',*args, **kwargs ):
	# 	'''	
	# 	anom = [numpy.ndarray] 2-d array representing a single monthly timestep of the data to be downscaled. 
	# 							Must also be representative of anomalies.
	# 	base = [str] filename of the corresponding baseline monthly file to use as template and downscale 
	# 							baseline for combining with anomalies.
	# 	src_transform = [affine.affine] 6 element affine transform of the input anomalies. [should be greenwich-centered]
	# 	resample_type = [str] one of ['bilinear', 'count', 'nearest', 'mode', 'cubic', 'index', 'average', 'lanczos', 'cubic_spline']
	# 	'''	
	# 	import rasterio
	# 	from rasterio.warp import reproject, RESAMPLING
	# 	from affine import Affine

	# 	resampling = {'average':RESAMPLING.average,
	# 				'cubic':RESAMPLING.cubic,
	# 				'lanczos':RESAMPLING.lanczos,
	# 				'bilinear':RESAMPLING.bilinear,
	# 				'cubic_spline':RESAMPLING.cubic_spline,
	# 				'mode':RESAMPLING.mode,
	# 				'count':RESAMPLING.count,
	# 				'index':RESAMPLING.index,
	# 				'nearest':RESAMPLING.nearest }
		
	# 	# # lets try to flip the data and affine and do this right.
	# 	# a,b,c,d,e,f,g,h,i = src_transform
	# 	# src_transform = Affine( a, b, c, d, -(e), np.abs(f) ) # DANGEROUS
	# 	# anom = np.flipud( anom )
	# 	# # end new stuff for flipping... <-- this should happen before the anoms and the src_transform get to this point.

	# 	base = rasterio.open( base )
	# 	baseline_arr = base.read( 1 )
	# 	baseline_meta = base.meta
	# 	baseline_meta.update( compress='lzw' )
	# 	output_arr = np.empty_like( baseline_arr )

	# 	reproject( anom, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata,
	# 				dst_transform=baseline_meta['affine'], dst_crs=baseline_meta['crs'],
	# 				dst_nodata=dst_nodata, resampling=resampling[ resample_type ], SOURCE_EXTRA=1000 )
		
	# 	return output_arr

# def sort_files( files, split_on='_', elem_month=-2, elem_year=-1 ):
# 	'''
# 	sort a list of files properly using the month and year parsed
# 	from the filename.  This is useful with SNAP data since the standard
# 	is to name files like '<prefix>_MM_YYYY.tif'.  If sorted using base
# 	Pythons sort/sorted functions, things will be sorted by the first char
# 	of the month, which makes thing go 1, 11, ... which sucks for timeseries
# 	this sorts it properly following SNAP standards as the default settings.
# 	ARGUMENTS:
# 	----------
# 	files = [list] list of `str` pathnames to be sorted by month and year. usually from glob.glob.
# 	split_on = [str] `str` character to split the filename on.  default:'_', SNAP standard.
# 	elem_month = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
# 		default:-2. For SNAP standard.
# 	elem_year = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
# 		default:-1. For SNAP standard.
# 	RETURNS:
# 	--------
# 	sorted `list` by month and year ascending. 
# 	'''
# 	import pandas as pd
# 	months = [ int(fn.split('.')[0].split( split_on )[elem_month]) for fn in files ]
# 	years = [ int(fn.split('.')[0].split( split_on )[elem_year]) for fn in files ]
# 	df = pd.DataFrame( {'fn':files, 'month':months, 'year':years} )
# 	df_sorted = df.sort_values( ['year', 'month' ] )
# 	return df_sorted.fn.tolist()

# def only_years( files, begin=1901, end=2100, split_on='_', elem_year=-1 ):
# 	'''
# 	return new list of filenames where they are truncated to begin:end
# 	ARGUMENTS:
# 	----------
# 	files = [list] list of `str` pathnames to be sorted by month and year. usually from glob.glob.
# 	begin = [int] four digit integer year of the begin time default:1901
# 	end = [int] four digit integer year of the end time default:2100
# 	split_on = [str] `str` character to split the filename on.  default:'_', SNAP standard.
# 	elem_year = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
# 		default:-1. For SNAP standard.
# 	RETURNS:
# 	--------
# 	sliced `list` to begin and end year.
# 	'''
# 	import pandas as pd
# 	years = [ int(fn.split('.')[0].split( split_on )[elem_year]) for fn in files ]
# 	df = pd.DataFrame( { 'fn':files, 'year':years } )
# 	df_slice = df[ (df.year >= begin ) & (df.year <= end ) ]
# 	return df_slice.fn.tolist()
