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
		self.filelist.sort()
		self.meta = rasterio.open( self.filelist[0] ).meta
		self.arrlist = [ rasterio.open( fn ).read( 1 ) for fn in self.filelist ]
		
		# make sure the data are sorted chronlogically
		# arrlist = self.arrlist
		# self.arrlist = self.sort_files()
		# self.rasters = [ rasterio.open(i).read(1) for i in self.arrlist ]

	# def sort_files( self ):
	# 	'''
	# 	ensure the files are in chronological order
	# 	'''
	# 	arrlist = self.arrlist
	# 	arrlist.sort()
	# 	self.arrlist = arrlist


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
					if False (default) do nothing.interp=False, ncpus=32,
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
		self._lonpc = None
		if interp:
			print( 'running interpolation across NAs' )
			_ = self.interp_na( )

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
	def interp_na( self ):
		'''
		np.float32
		method = [str] one of 'cubic', 'near', 'linear'

		return a list of dicts to pass to the xyz_to_grid hopefully in parallel
		'''
		from copy import copy
		# from pathos.multiprocessing import Pool
		# from multiprocessing import Pool
		from pathos.mp_map import mp_map
		import pandas as pd
		import numpy as np
		# remove the darn scientific notation
		np.set_printoptions( suppress=True )

		print 'interp with %s' % self.ncpus
		output_dtype = np.float32
		
		# if 0-360 leave it alone
		if ( self.ds.lon > 200.0 ).any() == True:

			dat, lons = self.ds[ self.variable ].data, self.ds.lon
			self._lonpc = lons
		else:
			# greenwich-centered rotate to 0-360 for interpolation across pacific
			dat, lons = self.rotate( self.ds[ self.variable ].data, self.ds.lon, to_pacific=True )
			self._rotated = True # update the rotated attribute
			self._lonpc = lons

		# mesh the lons and lats and unravel them to 1-D
		xi,yi = np.meshgrid( lons, self.ds.lat.data )
		lo, la = [ i.flatten() for i in (xi,yi) ]

		# setup args for multiprocessing
		df_list = [ pd.DataFrame({ 'x':lo, 'y':la, 'z':d.flatten() }).dropna( axis=0, how='any' ) for d in dat ]

		args = [ {'x':copy(np.array(df['x'])), 'y':copy(np.array(df['y'])), 'z':copy(np.array(df['z'])), \
				'grid':(copy(xi),copy(yi)), 'method':copy(self.method), 'output_dtype':copy(output_dtype) } for df in df_list ]

		# pool = Pool( self.ncpus )
		# out = pool.map( self._interpna, args[:50] )
		print( 'processing cru re-gridding in serial due to multiprocessing issues...' )
		out = map( self._interpna, args )
		# out = mp_map( self._interpna, args[:20], nproc=self.ncpus )
		# pool.close()
		lons = self._lonpc
		# stack em and roll-its axis so time is dim0
		dat = np.rollaxis( np.dstack( out ), -1 )
		if self._rotated == True: # rotate it back
			dat, lons = self.rotate( dat, lons, to_pacific=False )
				
		# place back into a new xarray.Dataset object for further processing
		ds = self.ds[ self.variable ] = dat
		var = ds[ self.variable ]
		setattr( var, 'data', dat )
		self.ds = ds
		print( 'ds interpolated updated into self.ds' )
		return 1
	@staticmethod
	def _interpna( args_dict ):
		return utils.xyz_to_grid( **args_dict )

