# # # # #
# Tool to downscale the CMIP5 data from the PCMDI group. 
# # # # #
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
			See Rasterio documentation for more details. 

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
	# make this a simple regrid command instead of interpolating the anomalies
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


class DownscaleAR5( object ):
	def __init__( self, ar5_modeled=None, ar5_historical=None, base_path=None, clim_path=None, climatology_begin='1961', climatology_end='1990', \
		plev=None, absolute=True, metric='metric', variable=None, ncores=2, *args, **kwargs ):
		'''
		NEW METHODS FOR AR5 DOWNSCALING USING THE NEW
		API-ECOSYSTEM.
		'''
		self.ar5_modeled = ar5_modeled
		self.ar5_historical = ar5_historical
		self.base_path = base_path
		self.clim_path = clim_path
		self.climatology_begin = climatology_begin
		self.climatology_end = climatology_end
		self.plev = plev
		self.absolute = absolute
		self.metric = metric
		self.variable = variable
		self.ncores = ncores
		self.utils = DownscalingUtils()

	@staticmethod
	def standardized_fn_to_vars( fn ):
		''' take a filename string following the convention for this downscaling and break into parts and return a dict'''
		name_convention = [ 'variable', 'cmor_table', 'model', 'scenario', 'experiment', 'begin_time', 'end_time' ]
		fn = os.path.basename( fn )
		fn_list = fn.split( '.' )[0].split( '_' )
		return { i:j for i,j in zip( name_convention, fn_list ) }
	def _calc_anomalies( self, *args, **kwargs ):
		'''
		calculate absolute or relative anomalies given a NetCDF file
		of the Climatic Research Unit (CRU) Historical Time Series.
		'''
		import xray

		# handle modeled vs. historical
		if self.ar5_modeled != None and self.ar5_historical != None:
			# parse the input name for some file metadata HARDWIRED!
			output_naming_dict = DownscaleAR5.standardized_fn_to_vars( self.ar5_modeled )
			variable = output_naming_dict[ 'variable' ]

			# read in both modeled and historical
			ds = xray.open_dataset( self.ar5_modeled )
			ds = ds[ variable ]
			clim_ds = xray.open_dataset( self.ar5_historical )
			# climatology
			clim_ds = clim_ds.loc[ {'time':slice(self.climatology_begin,self.climatology_end)} ]
			climatology = clim_ds[ variable ].groupby( 'time.month' ).mean( 'time' )
			del clim_ds

		elif self.ar5_historical is not None and self.ar5_modeled is None:
			output_naming_dict = standardized_fn_to_vars( self.ar5_historical )
			variable = output_naming_dict[ 'variable' ]

			# read in historical
			ds = xray.open_dataset( self.ar5_historical )
			# climatology
			climatology = ds.loc[ {'time':slice(self.climatology_begin,self.climatology_end)} ]
			climatology = climatology[ variable ].groupby( 'time.month' ).mean( 'time' )

		else:
			NameError( 'ERROR: must have both ar5_modeled and ar5_historical, or just ar5_historical' )

		if self.plev is not None:
			plevel, = np.where( ds.plev == self.plev )
			ds = ds[ :, plevel[0], ... ]
			climatology = climatology[ :, plevel[0], ... ]

		# anomalies
		if self.absolute == True:
			anomalies = ds.groupby( 'time.month' ) - climatology
		elif self.absolute == False:
			anomalies = ds.groupby( 'time.month' ) / climatology
		else:
			AttributeError( '_calc_anomalies (ar5): absolute can only be True or False' )
		return anomalies
	
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
		elif isinstance( anom, tuple ):
			anom_arr, meta = anom
		else:
			AttributeError( '_interp_downscale_wrapper: passed wrong instance type' )

		args_dict.update( output_filename=output_filename, anom_arr=anom_arr, meta=meta )
		return self.utils.downscale( **args_dict )

	def downscale_ar5_ts( self, *args, **kwargs ):
		from pathos.mp_map import mp_map

		# build output dirs

		# template setup

		# calc the anomalies
		anomalies = self._calc_anomalies()
		anomalies_pcll, lons_pcll = self.utils.shiftgrid( 0., anomalies, anomalies.lon.data ) # grabs lons from the xray ds

		# mesh the lons and lats and unravel them to 1-D
		lo, la = [ i.ravel() for i in np.meshgrid( lons_pcll, anomalies.lat ) ]
		
		# convert into pandas.DataFrame and drop all the NaNs -- land-only dataset
		anom_df_list = [ pd.DataFrame({ 'anom':i.ravel(), 'lat':la, 'lon':lo }).dropna( axis=0, how='any' ) for i in anomalies_pcll ]
		xi, yi = np.meshgrid( lons_pcll, anomalies.lat.data )

		# argument setup -- HARDWIRED
		# src_transform = affine.Affine( 0.5, 0.0, -180.0, 0.0, -0.5, 90.0 )
		# src_nodata = -9999.0
		# [!] THE ABOVE ARE INCORRECT FOR THE MODELED DATA


		# output_filenames setup
		dates = ds.time.to_pandas()
		years = dates.apply( lambda x: x.year ).tolist()
		months = [ i if len(i)==2 else '0'+i for i in np.arange( 1, 12+1, 1 ).astype( str ).tolist() ]
		month_year = [ (month, year) for year in years for month in months ]

		# read in the pre-processed 12-month climatology
		clim_list = sorted( glob.glob( os.path.join( self.clim_path, '*.tif' ) ) ) # this could catch you.
		clim_dict = { month:rasterio.open( fn ).read( 1 ) for month, fn in zip( months, clim_list ) }
		# [!] THIS BELOW NEEDS RE-WORKING FOR THE AR5 DATA MODELED DATA
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
		# ( anom_df, meshgrid_tuple, template_raster_fn, lons_pcll, src_transform, src_crs, src_nodata, output_filename, write_anomalies ) 	
		out = mp_map( lambda args: self._interp_downscale_wrapper( args_dict=args ), args_list, nproc=self.ncores )
		return 'downscaling complete. files output at: %s' % base_path

if __name__ == '__main__':
	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	# example of use of the new DownscaleCRU / DownscalingUtils classes
	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	import os, rasterio, xray, glob
	import pandas as pd
	import numpy as np

	# input args
	ar5_modeled = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped/IPSL-CM5A-LR/clt/clt_Amon_IPSL-CM5A-LR_rcp26_r1i1p1_200601_210012.nc'
	ar5_historical = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped/IPSL-CM5A-LR/clt/clt_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001_200512.nc'
	clim_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
	template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	base_path = '/atlas_scratch/malindgren/CMIP5'

	# EXAMPLE RUN -- TESTING
	down = DownscaleAR5( ar5_modeled, ar5_historical, base_path, clim_path, ncores=32) #, climatology_begin, climatology_end, plev, absolute, metric, ncores )
	output = down.downscale_ar5_ts()


# # # OLD SHIT BELOW 


# if __name__ == '__main__':
# 	import pandas as pd
# 	import numpy as np
# 	import os, sys, re, xray, rasterio, glob, argparse
# 	from rasterio import Affine as A
# 	from rasterio.warp import reproject, RESAMPLING
# 	from pathos import multiprocessing as mp

# 	# parse the commandline arguments
# 	parser = argparse.ArgumentParser( description='preprocess cmip5 input netcdf files to a common type and single files' )
# 	parser.add_argument( "-mi", "--modeled_fn", nargs='?', const=None, action='store', dest='modeled_fn', type=str, help="path to modeled input filename (NetCDF); default:None" )
# 	parser.add_argument( "-hi", "--historical_fn", nargs='?', const=None, action='store', dest='historical_fn', type=str, help="path to historical input filename (NetCDF); default:None" )
# 	parser.add_argument( "-o", "--output_path", action='store', dest='output_path', type=str, help="string path to the output folder containing the new downscaled outputs" )
# 	parser.add_argument( "-cbt", "--climatology_begin_time", nargs='?', const='196101', action='store', dest='climatology_begin', type=str, help="string in format YYYYMM or YYYY of the beginning month and potentially (year) of the climatology period" )
# 	parser.add_argument( "-cet", "--climatology_end_time", nargs='?', const='199012', action='store', dest='climatology_end', type=str, help="string in format YYYYMM or YYYY of the ending month and potentially (year) of the climatology period" )
# 	parser.add_argument( "-plev", "--plev", nargs='?', const=None, action='store', dest='plev', type=int, help="integer value (in millibars) of the desired pressure level to extract, if there is one." )
# 	parser.add_argument( "-cru", "--cru_path", action='store', dest='cru_path', type=str, help="path to the directory storing the cru climatology data derived from CL2.0" )
# 	parser.add_argument( "-at", "--anomalies_calc_type", nargs='?', const='absolute', action='store', dest='anomalies_calc_type', type=str, help="string of 'proportional' or 'absolute' to inform of anomalies calculation type to perform." )
# 	parser.add_argument( "-m", "--metric", nargs='?', const='metric', action='store', dest='metric', type=str, help="string of whatever the metric type is of the outputs to put in the filename." )
# 	parser.add_argument( "-dso", "--downscale_operation", action='store', dest='downscale_operation', type=str, help="string of 'add', 'mult', 'div', which refers to the type or downscaling operation to use." )
# 	parser.add_argument( "-nc", "--ncores", nargs='?', const=2, action='store', dest='ncores', type=int, help="integer valueof number of cores to use. default:2" )

# 	# parse args
# 	args = parser.parse_args()

# 	# unpack args
# 	modeled_fn = args.modeled_fn
# 	historical_fn = args.historical_fn
# 	output_path = args.output_path
# 	climatology_begin = args.climatology_begin
# 	climatology_end = args.climatology_end
# 	plev = args.plev
# 	cru_path = args.cru_path
# 	anomalies_calc_type = args.anomalies_calc_type
# 	metric = args.metric
# 	downscale_operation = args.downscale_operation
# 	ncores = args.ncores



# # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# 	# THIS APPEARS TO BE THE MAIN NOT THE DOWNSCALER.
# 	def downscale( src, dst, cru, src_crs, src_affine, dst_crs, dst_affine, output_filename, dst_meta, variable,\
# 		method='cubic_spline', operation='add', output_dtype='float32', **kwargs ):
# 		'''
# 		operation can be one of two keywords for the operation to perform the delta downscaling
# 		- keyword strings are one of: 'add'= addition, 'mult'=multiplication, or 'div'=division (not implemented)
# 		- method can be one of 'cubic_spline', 'nearest', 'bilinear' and must be input as a string.
# 		- output_dtype can be one of 'int32', 'float32'
# 		'''
# 		from rasterio.warp import reproject, RESAMPLING
# 		def add( cru, anom ):
# 			return cru + anom
# 		def mult( cru, anom ):
# 			return cru * anom
# 		def div( cru, anom ):
# 			# return cru / anom
# 			# this one may not be useful, but the placeholder is here 
# 			return NotImplementedError

# 		# switch to deal with numeric output dtypes
# 		dtypes_switch = {'int32':np.int32, 'float32':np.float32}

# 		# switch to deal with different resampling types
# 		method_switch = { 'nearest':RESAMPLING.nearest, 'bilinear':RESAMPLING.bilinear, 'cubic_spline':RESAMPLING.cubic_spline }
# 		method = method_switch[ method ]

# 		# reproject src to dst
# 		out = np.zeros( dst.shape ) 
# 		reproject( src,
# 					out,
# 					src_transform=src_affine,
# 					src_crs=src_crs,
# 					dst_transform=dst_affine,
# 					dst_crs=dst_crs,
# 					resampling=method )
# 		# switch to deal with different downscaling operators
# 		operation_switch = { 'add':add, 'mult':mult, 'div':div }
# 		downscaled = operation_switch[ operation ]( cru, out )

# 		# reset any > 100 values to 95 if the variable is cld or hur
# 		if variable == 'clt' or variable == 'hur' or variable == 'cld':
# 			downscaled[ downscaled > 100.0 ] = 95.0

# 		# give the proper fill values to the oob regions
# 		downscaled.fill_value = dst_meta['nodata']
# 		downscaled = downscaled.filled()

# 		# this is a geotiff creator so lets pass in the lzw compression
# 		dst_meta.update( compress='lzw' )
# 		with rasterio.open( output_filename, 'w', **dst_meta ) as out:
# 			out.write( downscaled.astype( dtypes_switch[ output_dtype ] ), 1 )
# 		return output_filename









	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	# [NOTE]: hardwired raster metadata meeting the ALFRESCO Model's needs for 
	# perfectly aligned inputs this is used as template metadata that 
	# is used in output generation. template raster filename below:
	# '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/
	#	TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	 # NO! 
	# meta_3338 = {'affine': A(2000.0, 0.0, -2173223.206087799, 
	# 				0.0, -2000.0, 2548412.932644147),
	# 			'count': 1,
	# 			'crs': {'init':'epsg:3338'},
	# 			'driver': u'GTiff',
	# 			'dtype': 'float32',
	# 			'height': 1186,
	# 			'nodata': -3.4e+38,
	# 			'width': 3218,
	# 			'compress':'lzw'}

	# # output template numpy array same dimensions as the template
	# dst = np.empty( (1186, 3218) )
	
	# # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	# # condition to deal with reading in historical data if needed.
	# if modeled_fn != None and historical_fn != None:
	# 	# parse the input name for some file metadata
	# 	output_naming_dict = standardized_fn_to_vars( modeled_fn )

	# 	# this is to maintain cleanliness
	# 	variable = output_naming_dict[ 'variable' ]

	# 	# read in both modeled and historical
	# 	ds = xray.open_dataset( modeled_fn )
	# 	ds = ds[ variable ].load()
	# 	clim_ds = xray.open_dataset( historical_fn )
	# 	clim_ds = clim_ds[ variable ].load()
	# 	# generate climatology / anomalies
	# 	clim_ds = clim_ds.loc[ {'time':slice(climatology_begin,climatology_end)} ]
	# 	climatology = clim_ds.groupby( 'time.month' ).mean( 'time' )

	# 	# find the begin/end years of the prepped files
	# 	dates = ds.time.to_pandas()
	# 	years = dates.apply( lambda x: x.year )
	# 	begin_time = years.min()
	# 	end_time = years.max()

	# 	del clim_ds
	# elif historical_fn is not None and modeled_fn is None:
	# 	# parse the input name for some file metadata
	# 	output_naming_dict = standardized_fn_to_vars( historical_fn )
		
	# 	# this is to maintain cleanliness
	# 	variable = output_naming_dict[ 'variable' ]

	# 	# read in historical
	# 	ds = xray.open_dataset( historical_fn )
	# 	ds = ds[ variable ].load()
	# 	# generate climatology / anomalies
	# 	climatology = ds.loc[ {'time':slice(climatology_begin,climatology_end)} ]
	# 	climatology = climatology.groupby( 'time.month' ).mean( 'time' )

	# 	# find the begin/end years of the prepped files
	# 	dates = ds.time.to_pandas()
	# 	years = dates.apply( lambda x: x.year )
	# 	begin_time = years.min()
	# 	end_time = years.max()

	# else:
	# 	NameError( 'ERROR: must have both modeled_fn and historical_fn, or just historical_fn' )

	# standardize the output pathing
	if output_naming_dict[ 'variable' ] == 'clt':
		variable_out = 'cld'
	else:
		variable_out = output_naming_dict[ 'variable' ]

	output_path = os.path.join( output_path, 'ar5', output_naming_dict['model'], variable_out, 'downscaled' )
	if not os.path.exists( output_path ):
		os.makedirs( output_path )


	# # if there is a pressure level to extract, extract it
	# if plev is not None:
	# 	plevel, = np.where( ds.plev == plev )
	# 	ds = ds[ :, plevel[0], ... ]
	# 	climatology = climatology[ :, plevel[0], ... ]

	# deal with different anomaly calculation types
	if anomalies_calc_type == 'absolute':
		anomalies = ds.groupby( 'time.month' ) - climatology
	elif anomalies_calc_type == 'proportional':
		anomalies = ds.groupby( 'time.month' ) / climatology
	else:
		NameError( 'anomalies_calc_type can only be one of "absolute" or "proportional"' )

	# some setup of the output raster metadata
	time_len, rows, cols = anomalies.shape
	crs = 'epsg:4326'
	affine = A( *[np.diff( ds.lon )[ 0 ], 0.0, -180.0, 0.0, -np.diff( ds.lat )[ 0 ], 90.0] )
	count = time_len
	resolution = ( np.diff( ds.lat )[ 0 ], np.diff( ds.lon )[ 0 ] )

	# close the dataset and clean it up
	ds = None

	# shift the grid to Greenwich Centering
	dat, lons = shiftgrid( 180., anomalies[:], anomalies.lon.data, start=False )

	# metadata for input?
	meta_4326 = {'affine':affine,
				'height':rows,
				'width':cols,
				'crs':crs,
				'driver':'GTiff',
				'dtype':np.float32,
				'count':time_len,
				'compress':'lzw' }
	# build some filenames for the outputs to be generated
	# months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
	months = [ i if len(i)==2 else '0'+i for i in np.arange( 1, 12+1, 1 ).astype( str ).tolist() ]

	years = [ str(year) for year in range( begin_time, end_time + 1, 1 ) ]
	# combine the months and the years
	combinations = [ (month, year) for year in years for month in months ]

	output_filenames = [ os.path.join( output_path, '_'.join([variable_out, 'metric', output_naming_dict['model'], output_naming_dict['scenario'], output_naming_dict['experiment'], month, year]) + '.tif' ) for month, year in combinations ]

	# load the baseline CRU CL2.0 data 
	# [NOTE]: THIS ASSUMES THEY ARE THE ONLY FILES IN THE DIRECTORY -- COULD BE A GOTCHA
	cru_files = glob.glob( os.path.join( cru_path, '*.tif' ) )
	cru_files.sort()
	cru_stack = [ rasterio.open( fn ).read( 1 ) for fn in cru_files ]
	# this is a hack to make a masked array with the cru data
	cru_stack = [ np.ma.masked_where( cru == cru.min(), cru ) for cru in cru_stack ]
	import itertools

	cru_gen = clim_generator( len(output_filenames), cru_stack )

	# cleanup some uneeded vars that are hogging RAM
	del climatology, anomalies

	# run in parallel using PATHOS
	pool = mp.Pool( processes=ncores )
	args_list = 

	args_list = [{ 'src':src, 
				'output_filename':fn, 
				'dst':dst, 
				'cru':cru, 
				'src_crs':meta_4326[ 'crs' ], 
				'src_affine':meta_4326[ 'affine' ],			
				'dst_crs':meta_3338[ 'crs' ], 
				'dst_affine':meta_3338[ 'affine' ], 
				'dst_meta':meta_3338, 
				'operation':downscale_operation, 
				'variable':variable }
				for src,fn,cru in zip( np.vsplit( dat, time_len ), output_filenames, cru_gen ) ]

	del dat, cru_gen, cru_stack

	out = pool.map( run, args_list )
	pool.close()
