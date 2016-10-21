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
		self.filelist = filelist
		self.meta = rasterio.open( self.filelist[0] ).meta
		self.arrlist = ( rasterio.open( fn ).read( 1 ) for fn in self.filelist )

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

		# update the lats and data to be NorthUp if necessary
		self._northup()

		# slice to the years we want if given
		if self.begin != None and self.end != None:
			self.ds = self.ds.sel( time=slice( str( self.begin ), str( self.end ) ) )

		self.interp = interp
		self.ncpus = ncpus
		self.method = method
		self.transform_from_latlon = utils.transform_from_latlon

	def _calc_affine( self ):
		return self.transform_from_latlon( self.ds.lat, self.ds.lon )
	def _northup( self, latitude='lat' ):
		''' this works only for global grids to be downscaled flips it northup '''
		if self.ds[ latitude ][0].data < 0: # meaning that south is north globally
			self.ds[ latitude ] = np.flipud( self.ds[ latitude ] )
			# flip each slice of the array and make a new one
			flipped = np.array( [ np.flipud( arr ) for arr in self.ds[ self.variable ].data ] )
			self.ds[ self.variable ] = (('time', 'lat', 'lon' ), flipped )
