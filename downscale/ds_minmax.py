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
import os
import numpy as np
import xarray as xr

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
		# setup new args
		self.mean_ds = mean_ds
		self.mean_variable = mean_variable

		super( DeltaDownscaleMinMax, self ).__init__( *args, **kwargs )
		
		# if there is no mean dataset to work with --> party's over
		if mean_ds == None:
			raise Exception( 'you must include the mean variable in the raw resolution \
								as arg `mean_ds`=downscale.Dataset object or use `DeltaDownscale`' )
	def _calc_climatolgy( self ):
		''' MASK THIS FOR MINMAX slice / aggregate to climatology using mean'''
		pass
	def _calc_anomalies( self ):
		''' calculate deltas but call them anomalies to fit the `downscale` pkg methods '''			
		self.anomalies = (self.historical.ds[ self.historical.variable ] - self.mean_ds.ds[ self.mean_variable ] ) #.to_dataset( name=variable )
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
		
		# rotate to greenwich-centered
		if ( self.anomalies.lon.data > 200.0 ).any() == True:
			dat, lons = utils.shiftgrid( 180., self.anomalies, self.anomalies.lon, start=False )
			self.anomalies_rot = dat
			a,b,c,d,e,f,g,h,i = self.affine
			# flip it to the greenwich-centering
			src_transform = self.historical.transform_from_latlon( self.historical.ds.lat, lons )
			# src_transform = affine.Affine( a, b, -180.0, d, e, 90.0 )
			print( 'anomalies rotated!' )
		else:
			dat, lons = ( self.anomalies, self.anomalies.lon )
			self.anomalies_rot = dat
			src_transform = self.historical.transform_from_latlon( self.historical.ds.lat, lons )
			# src_transform = Affine(0.5, 0.0, -180.0, 0.0, -0.5, 90.0)
			# src_transform = self.affine
			print( 'anomalies NOT rotated!' )

		# run and output
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
		f = partial( self.interp_ds, src_crs=self.src_crs, src_nodata=self.src_nodata, \
					dst_nodata=self.dst_nodata, src_transform=src_transform, resample_type=self.resample_type )

		run = partial( self._run_ds, f=f, operation_switch=operation_switch, anom=self.anom, mask_value=self.mask_value )

		# run it
		out = mp_map( run, args, nproc=self.ncpus )
		return output_dir

