# # #
# Downscale CRU Historical TS3.x data to a pre-processed climatology
#  extent, resolution, reference system
#
# Author: Michael Lindgren (malindgren@alaska.edu)
# # #

# import some modules
import rasterio, xray, os
import numpy as np
import pandas as pd
import numpy as np

class DownscalingUtils( object ):
	def write_gtiff( self, output_arr, template_meta, output_filename, compress=True ):
		'''
		DESCRIPTION:
		------------
		output a GeoTiff given a numpy ndarray, rasterio-style 
		metadata dictionary, and and output_filename.

		If a multiband file is to be processed, the Longitude
		dimension is expected to be the right-most. 
		--> dimensions should be (band, latitude, longitude)

		ARGUMENTS:
		----------
		output_arr = [numpy.ndarray] with longitude as the right-most dimension
		template_meta = [dict] rasterio-style raster meta dictionary.  Typically 
			found in a template raster by: rasterio.open( fn ).meta
		output_filename = [str] path to and name of the output GeoTiff to be 
			created.  currently only 'GTiff' is supported.
		compress = [bool] if True (default) LZW-compression is applied to the 
			output GeoTiff.  If False, no compression is applied.
			* this can also be added (along with many other gdal creation options)
			to the template meta as a key value pair template_meta.update( compress='lzw' ).
			See Rasterio documentation for more details. This is just a common one that is 

		RETURNS:
		--------
		string path to the new output_filename created

		'''
		import os
		if 'transform' in template_meta.keys():
			_ = template_meta.pop( 'transform' )
		if not output_filename.endswith( '.tif' ):
			UserWarning( 'output_filename does not end with ".tif", it has been fixed for you.' )
			output_filename = os.path.splitext( output_filename )[0] + '.tif'
		if output_arr.ndim == 2:
			# add in a new dimension - can get you into trouble with very large rasters...
			output_arr = output_arr[ np.newaxis, ... ] 
		elif output_arr.ndim < 2:
			raise ValueError( 'output_arr must have at least 2 dimensions' )
		nbands, nrows, ncols = output_arr.shape 
		if template_meta[ 'count' ] != nbands:
			raise ValueError( 'template_meta[ "count" ] must match output_arr bands' )
		if compress == True and 'compress' not in template_meta.keys():
			template_meta.update( compress='lzw' )
		with rasterio.open( output_filename, 'w', **template_meta ) as out:
			for band in range( 1, nbands+1 ):
				out.write( output_arr[ band-1, ... ], band )
		return output_filename
	def shiftgrid( self, lon0, datain, lonsin, start=True, cyclic=360.0 ):
		"""
		Shift global lat/lon grid east or west.
		.. tabularcolumns:: |l|L|
		==============   ====================================================
		Arguments        Description
		==============   ====================================================
		lon0             starting longitude for shifted grid
						 (ending longitude if start=False). lon0 must be on
						 input grid (within the range of lonsin).
		datain           original data with longitude the right-most
						 dimension.
		lonsin           original longitudes.
		==============   ====================================================
		.. tabularcolumns:: |l|L|
		==============   ====================================================
		Keywords         Description
		==============   ====================================================
		start            if True, lon0 represents the starting longitude
						 of the new grid. if False, lon0 is the ending
						 longitude. Default True.
		cyclic           width of periodic domain (default 360)
		==============   ====================================================
		returns ``dataout,lonsout`` (data and longitudes on shifted grid).
		"""
		if np.fabs(lonsin[-1]-lonsin[0]-cyclic) > 1.e-4:
			# Use all data instead of raise ValueError, 'cyclic point not included'
			start_idx = 0
		else:
			# If cyclic, remove the duplicate point
			start_idx = 1
		if lon0 < lonsin[0] or lon0 > lonsin[-1]:
			raise ValueError('lon0 outside of range of lonsin')
		i0 = np.argmin(np.fabs(lonsin-lon0))
		i0_shift = len(lonsin)-i0
		if np.ma.isMA(datain):
			dataout  = np.ma.zeros(datain.shape,datain.dtype)
		else:
			dataout  = np.zeros(datain.shape,datain.dtype)
		if np.ma.isMA(lonsin):
			lonsout = np.ma.zeros(lonsin.shape,lonsin.dtype)
		else:
			lonsout = np.zeros(lonsin.shape,lonsin.dtype)
		if start:
			lonsout[0:i0_shift] = lonsin[i0:]
		else:
			lonsout[0:i0_shift] = lonsin[i0:]-cyclic
		dataout[...,0:i0_shift] = datain[...,i0:]
		if start:
			lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]+cyclic
		else:
			lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]
		dataout[...,i0_shift:] = datain[...,start_idx:i0+start_idx]
		return dataout,lonsout
	def bounds_to_extent( self, bounds ):
		'''
		take input rasterio bounds object and return an extent
		'''
		l,b,r,t = bounds
		return [ (l,b), (r,b), (r,t), (l,t), (l,b) ]
	def padded_bounds( self, rst, npixels, crs ):
		'''
		convert the extents of 2 overlapping rasters to a shapefile with
		an expansion of the intersection of the rasters extents by npixels
		rst1: rasterio raster object
		rst2: rasterio raster object
		npixels: tuple of 4 (left(-),bottom(-),right(+),top(+)) number of pixels to
			expand in each direction. for 5 pixels in each direction it would look like
			this: (-5. -5. 5, 5) or just in the right and top directions like this:
			(0,0,5,5).
		crs: epsg code or proj4string defining the geospatial reference 
			system
		output_shapefile: string full path to the newly created output shapefile
		'''
		import rasterio, os, sys
		from shapely.geometry import Polygon

		resolution = rst.res[0]
		new_bounds = [ bound+(expand*resolution) for bound, expand in zip( rst.bounds, npixels ) ]
		return new_bounds
	def xyz_to_grid( self, x, y, z, grid, method='cubic', output_dtype=np.float32 ):
		'''
		interpolate points to a grid. simple wrapper around
		scipy.interpolate.griddata. Points and grid must be
		in the same coordinate system
		x = 1-D np.array of x coordinates / x,y,z must be same length
		y = 1-D np.array of y coordinates / x,y,z must be same length
		z = 1-D np.array of z coordinates / x,y,z must be same length
		grid = tuple of meshgrid as made using numpy.meshgrid()
				order (xi, yi)
		method = one of 'cubic', 'near', 'linear'
		'''
		from scipy.interpolate import griddata
		zi = griddata( (x, y), z, grid, method=method )
		zi = np.flipud( zi.astype( output_dtype ) )
		return zi
	def calc_anomalies( self, fn, variable, climatology_begin='1961', climatology_end='1990', absolute=True, *args, **kwargs ):
		'''
		calculate absolute or relative anomalies given a NetCDF file
		of the Climatic Research Unit (CRU) Historical Time Series.
		'''
		import xray
		ds = xray.open_dataset( fn )
		try:
			clim_ds = ds.loc[ {'time':slice(climatology_begin, climatology_end)} ]
			climatology = clim_ds[ variable ].groupby( 'time.month' ).mean( 'time' )
		except:
			AttributeError( 'cannot slice netcdf based on climatology years given. they must overlap.' )
		# calculate anomalies
		if absolute == True:
			anomalies = ds[ variable ].groupby( 'time.month' ) - climatology
		elif absolute == False:
			anomalies = ds[ variable ].groupby( 'time.month' ) / climatology
		else:
			AttributeError( 'calc_anomalies: absolute can only be True or False' )
		return anomalies
	def interpolate_anomalies( self, anom_df, meshgrid_tuple, template_raster_fn, lons_pcll, \
								src_transform, src_crs, src_nodata, output_filename, write_anomalies, *args, **kwargs ):
		'''
		run the interpolation to a grid, and reprojection / resampling to the Alaska / Canada rasters
		extent, resolution, origin (template_raster).

		This function is intended to be used to run a pathos.multiprocessing Pool's map function
		across a list of pre-computed arguments.

		ARGUMENTS:
		---------
		anom_df = []
		meshgrid_tuple = [] 
		template_raster_fn = [] 
		lons_pcll = [] 
		src_transform = [] 
		src_crs = [] 
		src_nodata = [] 
		output_filename = [] 
		write_anomalies = [] 
				
		RETURNS:
		-------

		if write_anomalies == True: [str] path to the output filename generated

		if write_anomalies == False: [tuple] interpolated NumPy ndarray representing the 
			interpolated anomalies and the rasterio-style metadata dictionary describing
			the newly generated raster.

		'''
		from rasterio.warp import reproject, RESAMPLING

		template_raster = rasterio.open( template_raster_fn )
		template_meta = template_raster.meta
		if 'transform' in template_meta.keys():
			template_meta.pop( 'transform' )
		# update some meta configs
		template_meta.update( compress='lzw', crs={'init':'epsg:3338'} )

		interp_arr = self.xyz_to_grid( np.array(anom_df['lon'].tolist()), \
						np.array(anom_df['lat'].tolist()), \
						np.array(anom_df['anom'].tolist()), grid=meshgrid_tuple, method='cubic' ) 

		src_nodata = -9999.0 # nodata
		interp_arr[ np.isnan( interp_arr ) ] = src_nodata
		dat, lons = self.shiftgrid( 180., interp_arr, lons_pcll, start=False )
		output_arr = np.empty_like( template_raster.read( 1 ) )

		reproject( dat, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
					dst_transform=template_meta['affine'], dst_crs=template_meta['crs'],\
					dst_nodata=None, resampling=RESAMPLING.cubic_spline, SOURCE_EXTRA=1000 )
		# mask it with the internal mask in the template raster, where 0 is oob.
		output_arr = np.ma.masked_where( template_raster.read_masks( 1 ) == 0, output_arr )
		output_arr.fill_value = template_meta[ 'nodata' ]
		output_arr = output_arr.filled()
		if write_anomalies == True:
			out = self.write_gtiff( output_arr, template_meta, output_filename, compress=True )
		elif write_anomalies == False:
			out = ( output_arr, template_meta )
		else:
			AttributeError( 'interpolate_anomalies: write_anomalies can be True or False only.')
		return out
	def downscale( self, anom_arr, baseline_arr, output_filename, \
					downscaling_operation, meta, post_downscale_function, *args, **kwargs ):
		'''
		downscale an anomaly array with a baseline array from the same period.

		Arguments:
		----------
		anom_arr = [ np.ndarray ] 2-D NumPy array representing a raster domain. 
					anom/baseline arrays must be same shape.
		baseline_arr = [ np.ndarray ] 2-D NumPy array representing a raster domain. 
					anom/baseline arrays must be same shape.
		output_filename = [ str ] full path and output filename to be created
		downscaling_operation = [ ] 
		meta = [ dict ] rasterio-style dictionary of raster metadata attributes. This 
				must jive with the dimensions and the data type of the array generated 
				through downscaling anom_arr with baseline_arr.  
		post_downscale_function = [ function ] a function that takes a 2-D downscaled 
				array as input and returns an array of the same shape / datatype.  This
				is typically used as a post-mortem for clamping the values from an output
				downscaled array that may be slightly outside the range due to the 
				interpolation method. We currently use this to clamp the values of the hur
				to 0-100.

		Returns:
		--------
		output_filename of newly generated downscaled raster.

		'''
		def add( base, anom ):
			return base + anom
		def mult( base, anom ):
			return base * anom
		def div( base, anom ):
			# this one may not be useful, but the placeholder is here
			# return base / anom
			return NotImplementedError

		try:
			operation_switch = { 'add':add, 'mult':mult, 'div':div }
		except:
			AttributeError( 'downscale: incorrect downscaling_operation str' )
		
		# [ CHECK ] This may be something better to be done before passing to this function
		# both files need to be masked here since we use a RIDICULOUS oob value...
		# for both tas and cld, values less than -200 are out of the range of acceptable values and it
		# grabs the -3.4... mask values. so lets mask using this
		baseline_arr = np.ma.masked_where( baseline_arr < -200, baseline_arr )
		anom_arr = np.ma.masked_where( anom_arr < -200, anom_arr )

		output_arr = operation_switch[ downscaling_operation ]( baseline_arr, anom_arr )
		output_arr[ np.isinf( output_arr ) ] = meta[ 'nodata' ]

		if post_downscale_function != None:
			output_arr = post_downscale_function( output_arr )

		if 'transform' in meta.keys():
			# avoid the gdal geotransform deprecation warning
			meta.pop( 'transform' )

		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( output_arr, 1 )
		return output_filename

class DownscaleCRU( object ):
	'''
	methods to downscale the Climatic Research Unit's (CRU) Historical 
	Time Series data using a 12-month climatology pre-processed to the final
	output domain and resolution.  Typically we use a PRISM climatology or a 
	CRU CL2.0 climatology for these purposes.

	'''
	def __init__( self, cru_ts, clim_path, template_raster_fn, base_path, climatology_begin='1961', climatology_end='1990', ncores=2, \
					absolute=True, metric='metric', variable=None, post_downscale_function=None, src_crs={'init':'epsg:4326'}, write_anomalies=True, *args, **kwargs ):
		self.cru_ts = cru_ts
		self.clim_path = clim_path
		self.template_raster_fn = template_raster_fn
		self.base_path = base_path
		self.climatology_begin = climatology_begin
		self.climatology_end = climatology_end
		self.ncores = ncores
		self.absolute = absolute
		self.metric = metric
		self.variable = variable
		self.post_downscale_function = post_downscale_function
		self.src_crs = src_crs
		self.utils = DownscalingUtils()
		self.write_anomalies = write_anomalies

	@staticmethod
	def _fn_month_grouper( fn, *args, **kwargs ):
		'''
		take a filename and return the month element of the naming convention
		'''
		return os.path.splitext( os.path.basename( fn ) )[0].split( '_' )[-2]
	def _get_varname_cru( self, *args, **kwargs ):
		'''
		take as input the cru ts3* netcdf filename and return (if possible)
		the name of the variable we want to work on from that netcdf.

		Arguments:
			nc_fn = [str] filepath to the cru ts* netcdf file used in downscaling

		Returns:
			the variable name as a string if it can be deduced, and errors if
			the variable name cannot be deduced.

		'''
		ds = xray.open_dataset( self.cru_ts )
		variables = ds.variables.keys()
		variable = [ variable for variable in variables \
						if variable not in [u'lon', u'lat', u'time'] ]
		if len( variable ) == 1:
			variable = variable[ 0 ]
		else:
			AttributeError( 'cannot deduce the variable from the file. supply nc_varname and re-run' )
		return variable
	def _get_years_cru( self, *args, **kwargs ):
		ds = xray.open_dataset( self.cru_ts )
		time = pd.DatetimeIndex( ds.time.values )
		years = [ year.year for year in time ]
		return years
	def _get_version_cru( self, *args, **kwargs ):
		version = ''.join( os.path.basename( self.cru_ts ).split( '.' )[:2] )
		version = version.replace( 'ts', 'TS' ) # to follow convention
		return version
	def _interp_downscale_wrapper( self, args_dict, *args, **kwargs  ):
		'''
		interpolate anomalies and downscale to the baseline arr
		'''
		output_filename = args_dict[ 'output_filename' ]
		args_dict.update( output_filename=output_filename.replace( 'downscaled', 'anom' ) )

		anom = self.utils.interpolate_anomalies( **args_dict )

		if isinstance( anom, basestring ):
			rst = rasterio.open( anom )
			meta = rst.meta
			meta.update( compress='lzw' )
			anom_arr = rst.read( 1 )
		elif isinstance( anom, tuple ): # isinstance( anom, tuple ):
			anom_arr, meta = anom
		else:
			AttributeError( '_interp_downscale_wrapper: passed wrong instance type' )

		args_dict.update( output_filename=output_filename, anom_arr=anom_arr, meta=meta )
		return self.utils.downscale( **args_dict )
	def downscale_cru_ts( self, *args, **kwargs ):
		'''
		run the CRU downscaling using the monthly climatology files given
		'''
		from pathos.mp_map import mp_map
		import glob, affine, rasterio

		nc_varname = self._get_varname_cru( )
		# handle cases where the desired varname != one parsed from file.
		if self.variable == None:
			variable = nc_varname
		else:
			variable = self.variable
		
		# build output dirs
		anomalies_path = os.path.join( base_path, variable, 'anom' )
		if not os.path.exists( anomalies_path ):
			os.makedirs( anomalies_path )

		downscaled_path = os.path.join( base_path, variable, 'downscaled' )
		if not os.path.exists( downscaled_path ):
			os.makedirs( downscaled_path )

		# template setup 
		template_raster = rasterio.open( self.template_raster_fn )
		template_meta = template_raster.meta
		template_meta.update( crs={'init':'epsg:3338'} )

		# make a mask with values of 0=nodata and 1=data
		template_raster_mask = template_raster.read_masks( 1 ) # mask of band 1 is all we need
		template_raster_mask[ template_raster_mask == 255 ] = 1

		anomalies = self.utils.calc_anomalies( self.cru_ts, variable, absolute=self.absolute )
		anomalies_pcll, lons_pcll = self.utils.shiftgrid( 0., anomalies, anomalies.lon.data ) # grabs lons from the xray ds

		# mesh the lons and lats and unravel them to 1-D
		lo, la = [ i.ravel() for i in np.meshgrid( lons_pcll, anomalies.lat ) ]
		
		# convert into pandas.DataFrame and drop all the NaNs -- land-only dataset
		anom_df_list = [ pd.DataFrame({ 'anom':i.ravel(), 'lat':la, 'lon':lo }).dropna( axis=0, how='any' ) for i in anomalies_pcll ]
		xi, yi = np.meshgrid( lons_pcll, anomalies.lat.data )

		# argument setup -- HARDWIRED
		src_transform = affine.Affine( 0.5, 0.0, -180.0, 0.0, -0.5, 90.0 )
		src_nodata = -9999.0
			
		# output_filenames setup
		years = np.unique( self._get_years_cru( self.cru_ts ) )
		cru_ts_version = self._get_version_cru( self.cru_ts ) # works if naming convention stays same
		months = [ i if len(i)==2 else '0'+i for i in np.arange( 1, 12+1, 1 ).astype( str ).tolist() ]
		month_year = [ (month, year) for year in years for month in months ]

		output_filenames = [ os.path.join( anomalies_path, '_'.join([ variable, self.metric, cru_ts_version, 'anom', month, str(year) ])+'.tif' )
								for month, year in month_year ]

		# NEW
		# read in the pre-processed 12-month climatology
		clim_list = sorted( glob.glob( os.path.join( self.clim_path, '*.tif' ) ) ) # this could catch you.
		clim_dict = { month:rasterio.open( fn ).read( 1 ) for month, fn in zip( months, clim_list ) }
		output_filenames = [ os.path.join( downscaled_path, '_'.join([ variable, self.metric, cru_ts_version, 'downscaled', month, str(year) ])+'.tif' )
								for month, year in month_year ]

		# set downscaling_operation based on self.absolute boolean
		if self.absolute == True:
			downscaling_operation = 'add'
		elif self.absolute == False:
			downscaling_operation = 'mult'
		else:
			AttributeError( 'downscaling operation: self.absolute must be boolean' )

		args_list = [ { 'anom_df':anom_df, 
						'meshgrid_tuple':(xi, yi), 
						'template_raster_fn':template_raster_fn, 
						'lons_pcll':lons_pcll, 
						'src_transform':src_transform, 
						'src_crs':self.src_crs, \
						'src_nodata':src_nodata, 
						'output_filename':out_fn,
						'baseline_arr':clim_dict[ self._fn_month_grouper( out_fn ) ],
						'downscaling_operation':downscaling_operation, 
						'post_downscale_function':self.post_downscale_function,
						'write_anomalies':self.write_anomalies }
							for anom_df, out_fn in zip( anom_df_list, output_filenames ) ]

		# run anomalies interpolation and downscaling in a single go.
		out = mp_map( lambda args: self._interp_downscale_wrapper( args_dict=args ), args_list, nproc=self.ncores )
		return 'downscaling complete. files output at: %s' % base_path

if __name__ == '__main__':
	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	# example of use of the new DownscaleCRU / DownscalingUtils classes
	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	# example of post_downscale_function - pass in at DownscaleCRU()
	def clamp_vals( x ):
		''' clamp the values following the relative humidity downscaling '''
		x[ (x > 100) & (x < 500) ] = 95
		return x

	# minimum required arguments
	cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.cld.dat.nc'
	clim_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
	template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	base_path = '/atlas_scratch/malindgren/CMIP5'
	
	# run example
	down = DownscaleCRU( cru_ts, clim_path, template_raster_fn, base_path, absolute=False, ncores=32 )
	output = down.downscale_cru_ts()

