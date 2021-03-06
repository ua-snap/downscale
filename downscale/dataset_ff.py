# xarray removed version -- due to some SEVERE issues in pd.TimeStamp
import rasterio, os
import numpy as np
import pandas as pd
import geopandas as gpd
from downscale import utils

class DatasetFF( object ):
	''' 
	THIS SHOULD BE SUBCLASSED FROM `Dataset`, but is not because that class
	is poorly written and the __init__ to complex (currently) to override without
	essentially rewriting the entire class (!TODO!).

	'''
	def __init__( self, fn, variable, model, scenario, project=None, units=None, metric=None, 
					interp=False, ncpus=32,	method='linear', level=None, level_name=None, begin=None, end=None, *args, **kwargs ):
		'''
		build a dataset object to store the NetCDF low-res data for downscaling.
		
		ARGUMENTS:
		----------
		fn = [str] str path name or list of str path names of netCDF4 dataset(s) to be read in -- RAW CMIP5 SORTED CHRONOLOGICALLY for MFDataset.
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
		level = [ NOT YET IMPLEMENTED ]
		level_name = [ NOT YET IMPLEMENTED ]
		begin = [int] desired begin year
		end = [int] desired end year

		'''
		import ast
		import netCDF4
		from netCDF4 import MFDataset, num2date

		self.fn = fn # CHRONOLOGICALLY SORTED! [list] of MFDataset-able filenames
		ds = MFDataset( self.fn )
		self.ds = ds
		try:
			t = ds.variables['time']
			dates = num2date( t[:], calendar=t.calendar, units=t.units )
			self.fileyear_begin = dates.min().year
			self.fileyear_end = dates.max().year
		except:
			raise AttributeError( '[downscale]: input netcdf datasets do not conform to CF-standards for \n\
							time, calendar, and/or (time) units' )

		self.variable = variable
		self.model = model
		self.scenario = scenario
		self.level = level
		self.level_name = level_name
		self.begin = begin
		self.end = end

		self.lon = ds['lon'][:]
		self.lat = ds['lat'][:]
		
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
		
		if self.begin is not None and self.end is not None:
			
			# for slicing to desired times and not what is in the files
			years = np.repeat( range(self.fileyear_begin, self.fileyear_end+1), 12)
			# begin_idx, = np.where( years == self.begin ) # ASSUMES FULL YEARS!!!!
			# begin_idx = int(begin_idx.min()) # months
			# end_idx, = np.where( years == self.end ) # ASSUMES FULL YEARS!!!!
			# end_idx = int(end_idx.max())+1 # months

			# will grab the closest years possible
			begin_idx = (np.abs(years - begin)).argmin()
			end_idx = (np.abs(years - end)).argmin() + 12 # months since will grab first instance.

			if self.level is not None and self.level_name is not None:
				NotImplementedError( 'LEVELS ARE NOT YET SUPPORTED BY FAR-FUTURES in `downscale`' )
				# THIS IS DANGEROUS! BUT WOULD BE HARD TO EXPLOIT DUE TO THE STRUCTURE
				# OF THE QUERY BEING WITHIN THIS MODULE...  Famous last words...
				ds_lev = eval( 'ds.sel({}={})'.format( self.level_name, self.level ) )
				# levidx, = np.where( ds[ self.level_name ] == self.level )
				# ds = ds[ self.variable ][ :, int(levidx), ... ]
				self.dat = ds_lev[self.variable][begin_idx:end_idx, ...]
				del ds_lev
			else:
				self.dat = ds[self.variable][begin_idx:end_idx, ...]
		else:
			self.dat = ds[self.variable][:]
			self.begin = self.fileyear_begin
			self.end = self.fileyear_end
		
		# update the lats and data to be NorthUp if necessary
		self._northup()

		self.interp = interp
		self.ncpus = ncpus
		self.method = 'linear'
		self.transform_from_latlon = utils.transform_from_latlon

	def _calc_affine( self ):
		''' 
		calculate affine transform from lats / lons and snap to global extent
		NOTE: only use for global data
		'''
		return self.transform_from_latlon( self.lat, self.lon )
	def _northup( self, latitude='lat' ):
		''' this works only for global grids to be downscaled flips it northup '''
		if self.lat[0] < 0: # meaning that south is north globally
			self.lat = np.flipud( self.lat )
			# flip each slice of the array and make a new one
			self.dat = np.array( [ np.flipud( arr ) for arr in self.dat ] )
			print( 'flipped to North-up' )


# # # keep temporarily
# # how to make a range of datetime64 obj's with numpy...
# climdates = np.arange('1961-01-01','1991-01-01', dtype='datetime64[M]')


class MaskFF( object ):
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
	def mask( self, latitude='lat', longitude='lon', all_touched=True ):
		''' make a mask from the aoi shapefile and the low-res input NetCDF '''
		import geopandas as gpd
		from downscale import utils

		gdf = gpd.read_file( self.aoi )
		shapes = [ (geom, self.mask_value) for geom in gdf.geometry ]
		ds = self.ds.ds # grab the ds sub-object from the Dataset object
		# coords = ds.coords # get lats and lons as a coords dict from xarray
		coords = self._get_coords() #[FF added] a hack to overcome our masking procedure when we can use xarray
		return utils.rasterize( shapes, coords=coords, latitude=latitude, longitude=longitude, fill=self.fill_value, **{'all_touched':all_touched} ).data
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
	def _get_coords( self ):
		import xarray as xr
		ds = xr.open_dataset(self.ds.fn[0], decode_times=False)
		coords = {'lon':np.array(ds.lon),'lat':np.flipud(ds.lat)}
		return coords

