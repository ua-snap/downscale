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
		self.arrlist = [ rasterio.open( fn ).read( 1 ) for fn in self.filelist ]


class Dataset( object ):
	def __init__( self, fn, variable, model, scenario, units=None, interp=False, ncpus=32, \
					method='cubic', *args, **kwargs):
		'''
		fn = [str] path to the xray supported dataset to be read in.
		variable = [str] abbreviation of variable name to extract from file
		model = [str] name of the model being read
		scenario = [str] name of the scenario being read
		units = [str] abbreviation of the units of the variable
		interp = [bool] if True interpolate across NA's using a spline. 
					if False (default) do nothing.
		ncpus = [ int ] number of cores to use if interp=True. default:2.
		'''
		import xarray as xr
		self.fn = fn
		self.ds = xr.open_dataset( fn )
		self.variable = variable
		self.model = model
		self.scenario = scenario
		self.units = units
		self.interp = interp
		self.ncpus = ncpus
		self.method = method
		self._rotated = False

		if interp:
			print( 'running interpolation across NAs' )
			# self._interpna( method=method )
			self.run = self.run()

	@staticmethod
	def rotate( dat, lons, to_pacific=False ):
		'''rotate longitudes in WGS84 Global Extent'''
		if to_pacific == False:
			# to -180.0 - 180.0 
			dat, lons = utils.shiftgrid( 180., dat, lons, start=False )
		elif to_pacific == True:
			# to 0 - 360
			dat, lons = utils.shiftgrid( 0., dat, lons )
		else:
			raise AttributeError( 'to_pacific can be only one of True or False' )
		return dat, lons
	@staticmethod
	def _wrap( x ):
		# a function to make this easier for function passing
		return utils.xyz_to_grid( **x )
	def _interpna_setup( self ):
		'''
		np.float32
		method = [str] one of 'cubic', 'near', 'linear'
		'''
		# import pathos
		# from pathos import multiprocessing
		# from pathos.mp_map import mp_map
		# from pathos.multiprocessing import ProcessingPool as Pool
		from functools import partial
		from copy import copy
		import multiprocessing
		print 'interp with %s' % self.ncpus
		output_dtype = np.float32
		
		# if greenwich-centered, lets rotate it to pcll since we are interested in ALASKA
		if ( self.ds.lon > 200.0 ).any() == True:
			dat, lons = self.ds[ self.variable ].data, self.ds.lon
			self._lonpc = lons
		else:
			dat, lons = self.rotate( self.ds[ self.variable ].data, self.ds.lon, to_pacific=True )
			self._rotated = True # update the rotated attribute
			self._lonpc = lons

		# mesh the lons and lats and unravel them to 1-D
		xi,yi = np.meshgrid( lons, self.ds.lat.data )
		lo, la = [ i.ravel() for i in (xi,yi) ]

		# setup args for multiprocessing
		args = [ pd.DataFrame( \
					{'x':copy(lo),'y':copy(la),'z':d.copy().ravel()} ).dropna( axis=0, how='any' ).to_dict( orient='list' ) \
					for d in dat ]
		_ = [ arg.update( grid=copy((xi,yi)), method=self.method, output_dtype=copy(output_dtype) ) for arg in args ]
		return args
	@staticmethod
	def _interpna( args_dict ):
		return utils.xyz_to_grid( **args_dict )
	def run( self ):
		# from pathos.multiprocessing import Pool
		from pathos.multiprocessing import ProcessPool as Pool
		args = self._interpna_setup( )
		pool = Pool( processes=self.ncpus )
		out = pool.map( self._interpna, args[:400] )
		pool.close()
		lons = self._lonpc
		# stack em and roll-its axis so time is dim0
		dat = np.rollaxis( np.dstack( out ), -1 )
		if self._rotated == True: # rotate it back
			dat, lons = self.rotate( dat, lons, to_pacific=False )
		# place back into a new xarray.Dataset object for further processing
		# function to make a new xarray.Dataset object with the mdata we need?
		# ds = self.ds
		# var = ds[ self.variable ]
		# setattr( var, 'data', dat )
		# self.ds = ds
		print( 'ds interpolated updated into self.ds' )
		return dat
