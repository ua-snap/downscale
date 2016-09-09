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
	def __init__( self, baseline, clim_begin, clim_end, historical, future, downscaling_operation, 
				level, level_name, mask, mask_value, ncpus, src_crs, src_nodata, dst_nodata, post_downscale_function, 
				varname, modelname, anom, resample_type, mean_ds=None, mean_variable=None ): # , *args, **kwargs
		'''
			note here that all data falls into the 'historical' category, because we no longer need to 
			have the 1961-1990 climatology period for the futures as this version of DeltaDownscale computes
			deltas by removing the mean in time instead of removing the climatology.
		'''
		print(1)
		super( DeltaDownscaleMinMax, self ).__init__( baseline, clim_begin, clim_end, historical, future, downscaling_operation, 
							level, level_name, mask, mask_value, ncpus, src_crs, src_nodata, dst_nodata, post_downscale_function, 
							varname, modelname, anom, resample_type )

		# setup new args
		self.mean_ds = mean_ds
		self.mean_variable = mean_variable

		# , mean_ds=None, mean_variable=None
		
		# if there is no mean dataset to work with --> party's over
		if mean_ds == None:
			raise Exception( 'you must include the mean variable in the raw resolution \
								as arg `mean_ds`=downscale.Dataset object or use `DeltaDownscale`' )

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



# # FOR RUN OF THE MIN / MAX TAS DATA:
# # 1. COMPUTE DELTAS FIRST ANND WRITE TO NETCDF
# # 2. USE `DeltaDownscaleMinMax` for the downscaling


# fn = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.tmx.dat.nc'
# mean_fn = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.tmp.dat.nc'
# variable = 'tmx'
# mean_variable = 'tmp'
# output_filename = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/test/cru_ts3.23.1901.2014.tmx_delta_tmp.dat.nc'
# _ = delta_mm( fn, mean_fn, variable, mean_variable, output_filename )

# # now use the new DeltaDownscaleMM class to do the work.



