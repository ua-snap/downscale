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
import geopandas as gpd
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

class Mask( object ):
	def __init__( self, aoi, ds, mask_value=1, fill_value=0, *args, **kwargs ):
		'''
		make a mask from a shapefile which is already in the CRS and domain of the 
		input ds.

		ARGUMENTS:
		---------
		aoi = [str] full read path of a shapefile with .shp extension
		ds = [downscale.Dataset] instance of file in a downscale.Dataset object
		mask_value = [int] value to use for masked areas. default:1.
		fill_value = [int] value to use for unmasked areas. default:0.

		'''
		self.aoi = aoi
		self.ds = ds
		self.mask_value = mask_value
		self.fill_value = fill_value
	@property
	def mask( self, latitude='lat', longitude='lon' ):
		''' make a mask from the aoi shapefile and the low-res input NetCDF '''
		import geopandas as gpd
		from downscale import utils

		gdf = gpd.read_file( self.aoi )
		shapes = [ (geom, self.mask_value) for geom in gdf.geometry ]
		ds = self.ds.ds # grab the ds sub-object from the Dataset object
		coords = ds.coords # get lats and lons as a coords dict from xarray
		return utils.rasterize( shapes, coords=coords, latitude=latitude, longitude=longitude, fill=self.fill_value ).data
	def to_gtiff( self, output_filename ):
		''' write the mask to geotiff given an output_filename '''
		meta = {'compress':'lzw'}
		count, height, width = self.ds.ds[ self.ds.variable ].shape
		affine = self.ds._calc_affine()
		meta.update( affine=affine, height=height, width=width, count=1, dtype='int16', driver='GTiff' )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( self.mask.astype( np.int16 ), 1 )
	def _dump_out_testfile( self, output_filename ):
		''' hidden function to dump out a single representative raster from the NetCDF for comparison '''
		meta = {'compress':'lzw'}
		count, height, width = self.ds.ds[ self.ds.variable ].shape
		affine = self.ds._calc_affine()
		meta.update( affine=affine, height=height, width=width, count=1, dtype='float', driver='GTiff' )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( self.ds.ds[self.ds.variable][0], 1 )

class Dataset( object ):
	def __init__( self, fn, variable, model, scenario, project=None, units=None, metric=None, 
					interp=False, ncpus=32,	method='linear', begin=None, end=None, *args, **kwargs ):
		'''
		build a dataset object to store the NetCDF low-res data for downscaling.
		
		ARGUMENTS:
		----------
		fn = [str] path to the xray supported dataset to be read in.
		variable = [str] abbreviation of variable name to extract from file
		model = [str] name of the model being read
		scenario = [str] name of the scenario being read
		project = [str] name of the project.  ex. 'ar5'
		units = [str] abbreviation of the units of the variable
		metric = [str] metric used to describe variable temporally. ex. 'mean', 'total'
		interp = [bool] if True interpolate across NA's using a spline. 
					if False (default) do nothing.interp=False, ncpus=32,
		ncpus = [ int ] number of cores to use if interp=True. default:2.
		method = [ str ] type of interpolation to use. hardwired to 'linear' currently
		begin = year series begins
		end = year series ends

		'''
		import xarray as xr
		self.fn = fn
		self.ds = xr.open_dataset( self.fn )
		self.variable = variable
		self.model = model
		self.scenario = scenario
		self.begin = begin # year begin
		self.end = end # year end
		
		if units:
			self.units = units
		else:
			self.units = 'units'
			
		if project:
			self.project = project
		else:
			self.project = 'project'

		if metric:
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
		self.method = 'linear'
		self.transform_from_latlon = utils.transform_from_latlon

	def _calc_affine( self ):
		''' 
		calculate affine transform from lats / lons and snap to global extent
		NOTE: only use for global data
		'''
		return self.transform_from_latlon( self.ds.lat, self.ds.lon )
	def _northup( self, latitude='lat' ):
		''' this works only for global grids to be downscaled flips it northup '''
		if self.ds[ latitude ][0].data < 0: # meaning that south is north globally
			self.ds[ latitude ] = np.flipud( self.ds[ latitude ] )
			# flip each slice of the array and make a new one
			flipped = np.array( [ np.flipud( arr ) for arr in self.ds[ self.variable ].data ] )
			self.ds[ self.variable ] = (('time', 'lat', 'lon' ), flipped )
