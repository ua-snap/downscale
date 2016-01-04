# # #
# Current implementation of the cru ts31 (ts32) delta downscaling procedure
#
# Author: Michael Lindgren (malindgren@alaska.edu)
# # #

# import some modules
import rasterio, xray, os
import numpy as np
import pandas as pd

class DownscaleCRU( object ):
	'''
	methods to downscale the Climatic Research Unit's (CRU) Historical 
	Time Series data using a 12-month climatology pre-processed to the final
	output domain and resolution.  Typically we use a PRISM climatology or a 
	CRU CL2.0 climatology for these purposes.

	'''
	def write_gtiff( output_arr, template_meta, output_filename, compress=True ):
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
	def shiftgrid( lon0, datain, lonsin, start=True, cyclic=360.0 ):
		import numpy as np
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
	def bounds_to_extent( bounds ):
		'''
		take input rasterio bounds object and return an extent
		'''
		l,b,r,t = bounds
		return [ (l,b), (r,b), (r,t), (l,t), (l,b) ]
	def padded_bounds( rst, npixels, crs ):
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
	def xyz_to_grid( x, y, z, grid, method='cubic', output_dtype=np.float32 ):
		'''
		interpolate points to a grid. simple wrapper around
		scipy.interpolate.griddata. Points and grid must be
		in the same coordinate system
		x = 1-D np.array of x coordinates / x,y,z must be same length
		y = 1-D np.array of y coordinates / x,y,z must be same length
		z = 1-D np.array of z coordinates / x,y,z must be same length
		grid = tuple of meshgrid as made using numpy.meshgrid()
				order (xi, yi)
		method = one of 'cubic', 'near', linear
		'''
		from scipy.interpolate import griddata
		zi = griddata( (x, y), z, grid, method=method )
		zi = np.flipud( zi.astype( output_dtype ) )
		return zi
	def _interpolate_anomalies( anomalies, meshgrid_tuple, lons_pcll, template_raster_fn, src_transform, src_crs, src_nodata, output_filename, *args, **kwargs ):
		'''
		run the interpolation to a grid, and reprojection / resampling to the Alaska / Canada rasters
		extent, resolution, origin (template_raster).

		This function is intended to be used to run a pathos.multiprocessing Pool's map function
		across a list of pre-computed arguments.
				
		RETURNS:

		[str] path to the output filename generated

		'''
		from rasterio.warp import reproject, RESAMPLING

		template_raster = rasterio.open( template_raster_fn )
		template_meta = template_raster.meta
		if 'transform' in template_meta.keys():
			template_meta.pop( 'transform' )
		# update some meta configs
		template_meta.update( crs={'init':'epsg:3338'} )
		template_meta.update( compress='lzw' )

		interp_arr = xyz_to_grid( np.array(anomalies['lon'].tolist()), \
						np.array(anomalies['lat'].tolist()), \
						np.array(anomalies['anom'].tolist()), grid=meshgrid_tuple, method='cubic' ) 

		src_nodata = -9999.0 # nodata
		interp_arr[ np.isnan( interp_arr ) ] = src_nodata
		dat, lons = shiftgrid( 180., interp_arr, lons_pcll, start=False )
		output_arr = np.empty_like( template_raster.read( 1 ) )

		reproject( dat, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
					dst_transform=template_meta['affine'], dst_crs=template_meta['crs'],\
					dst_nodata=None, resampling=RESAMPLING.cubic_spline, num_threads=1, SOURCE_EXTRA=1000 )	
		# mask it with the internal mask in the template raster, where 0 is oob.
		output_arr = np.ma.masked_where( template_raster.read_masks( 1 ) == 0, output_arr )
		output_arr.fill_value = template_meta[ 'nodata' ]
		output_arr = output_arr.filled()
		return write_gtiff( output_arr, template_meta, output_filename, compress=True )
	def _fn_month_grouper( x ):
		'''
		take a filename and return the month element of the naming convention
		'''
		return os.path.splitext(os.path.basename(x))[0].split( '_' )[-2]
	def _downscale_cru_historical( file_list, cru_cl20_arr, output_path, downscaling_operation ):
		'''
		take a list of cru_historical anomalies filenames, groupby month,
		then downscale with the cru_cl20 climatology as a numpy 2d ndarray
		that is also on the same grid as the anomalies files.
		(intended to be the akcan 1km/2km extent).

		operation can be one of 'mult', 'add', 'div' and represents the
		downscaling operation to be use to scale the anomalies on top of the baseline.
		this is based on how the anomalies were initially calculated.

		RETURNS:

		output path location of the new downscaled files.
		'''
		from functools import partial

		def f( anomaly_fn, baseline_arr, output_path, downscaling_operation ):
			def add( cru, anom ):
				return cru + anom
			def mult( cru, anom ):
				return cru * anom
			def div( cru, anom ):
				# return cru / anom
				# this one may not be useful, but the placeholder is here 
				return NotImplementedError

			cru_ts31 = rasterio.open( anomaly_fn )
			meta = cru_ts31.meta
			meta.update( compress='lzw', crs={'init':'epsg:3338'} )
			cru_ts31 = cru_ts31.read( 1 )
			operation_switch = { 'add':add, 'mult':mult, 'div':div }
			# this is hardwired stuff for this fairly hardwired script.
			output_filename = os.path.basename( anomaly_fn ).replace( 'anom', 'downscaled' )
			output_filename = os.path.join( output_path, output_filename )
			# both files need to be masked here since we use a RIDICULOUS oob value...
			# for both tas and cld, values less than -200 are out of the range of acceptable values and it
			# grabs the -3.4... mask values. so lets mask using this
			baseline_arr = np.ma.masked_where( baseline_arr < -200, baseline_arr )
			cru_ts31 = np.ma.masked_where( cru_ts31 < -200, cru_ts31 )

			output_arr = operation_switch[ downscaling_operation ]( baseline_arr, cru_ts31 )
			output_arr[ np.isinf(output_arr) ] = meta[ 'nodata' ]

			if post_downscale_function != None:
				output_arr = post_downscale_function( output_arr )

			if 'transform' in meta.keys():
				meta.pop( 'transform' )
			with rasterio.open( output_filename, 'w', **meta ) as out:
				out.write( output_arr, 1 )
			return output_filename
			
		partial_f = partial( f, baseline_arr=cru_cl20_arr, output_path=output_path, downscaling_operation=downscaling_operation )
		cru_ts31 = file_list.apply( lambda fn: partial_f( anomaly_fn=fn ) )
		return output_path
	# # # # # 
	def _get_varname_cru( nc_fn ):
		'''
		take as input the cru ts3* netcdf filename and return (if possible)
		the name of the variable we want to work on from that netcdf.

		Arguments:
			nc_fn = [str] filepath to the cru ts* netcdf file used in downscaling

		Returns:
			the variable name as a string if it can be deduced, and errors if
			the variable name cannot be deduced.

		'''
		ds = xray.open_dataset( nc_fn )
		variables = ds.variables.keys()
		variable = [ variable for variable in variables \
						if variable not in [u'lon', u'lat', u'time'] ]
		if len( variable ) == 1:
			variable = variable[ 0 ]
		else:
			AttributeError( 'cannot deduce the variable from the file. supply nc_varname and re-run' )
		return variable
	def _calc_anomalies( nc_fn, nc_varname, climatology_begin, climatology_end, absolute=True ):
		'''
		calculate absolute or relative anomalies given a NetCDF file
		of the Climatic Research Unit (CRU) Historical Time Series.
		'''
		ds = xray.open_dataset( nc_fn )
		# climatology -- slice the time dimension
		try:
			clim_ds = ds.loc[ {'time':slice(climatology_begin,climatology_end)} ]
			climatology = clim_ds[ nc_varname ].groupby( 'time.month' ).mean( 'time' )
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
	def _get_years_cru( nc_fn ):
		ds = xray.open_dataset( nc_fn )
		time = pd.DatetimeIndex( ds.time.values )
		years = [ year.year for year in time ]
		return years
	def _get_version_cru( nc_fn ):
		version = ''.join(os.path.basename( nc_fn ).split( '.' )[:2])
		version = version.replace( 'ts', 'TS' ) # to follow convention
		return version
	def _main( x, *args, **kwargs ):
		'''
		run the CRU downscaling using the monthly climatology files given
		'''
		from pathos.mp_map import mp_map
		import glob, affine

		nc_varname = get_varname_cru( nc_fn )
		# handle cases where the desired varname is not the same as the one parsed from file.
		if variable == None:
			variable = nc_varname
		else:
			variable = nc_varname
		
		# build output dirs
		anomalies_path = os.path.join( base_path, variable, 'anom' )
		if not os.path.exists( anomalies_path ):
			os.makedirs( anomalies_path )

		downscaled_path = os.path.join( base_path, variable, 'downscaled' )
		if not os.path.exists( downscaled_path ):
			os.makedirs( downscaled_path )

		# template setup 
		template_raster = rasterio.open( template_raster_fn )
		template_meta = template_raster.meta
		template_meta.update( crs={'init':'epsg:3338'} )

		# make a mask with values of 0=nodata and 1=data
		template_raster_mask = template_raster.read_masks( 1 ) # mask of band 1 is all we need
		template_raster_mask[ template_raster_mask == 255 ] = 1

		anomalies = calc_anomalies( nc_fn, nc_varname, climatology_begin, climatology_end, absolute ) # the absolute calculation needs some thought
		anomalies_pcll, lons_pcll = shiftgrid( 0., anomalies, anomalies.lon.data ) # grabs lons from the xray ds

		# mesh the lons and lats and unravel them to 1-D
		lo, la = [ i.ravel() for i in np.meshgrid( lons_pcll, anomalies.lat ) ]
		
		# convert into pandas.DataFrame and drop all the NaNs -- land-only dataset
		anom_df_list = [ pd.DataFrame({ 'anom':i.ravel(), 'lat':la, 'lon':lo }).dropna( axis=0, how='any' ) for i in dat_pcll ]
		xi, yi = np.meshgrid( lons_pcll, anomalies.lat.data )

		# argumet setup -- HARDWIRED
		src_transform = affine.Affine( 0.5, 0.0, -180.0, 0.0, -0.5, 90.0 )
		src_crs = {'init':'epsg:4326'}
		src_nodata = -9999.0
			
		# output_filenames setup
		years = get_years_cru( nc_fn )
		cru_ts_version = get_version_cru( nc_fn ) # works if naming convention stays same
		months = [ i if len(i)==2 else '0'+i for i in np.arange( 1, 12+1, 1 ).astype( str ).tolist() ]
		month_year = [ (month, year) for year in years for month in months ]

		output_filenames = [ os.path.join( anomalies_path, '_'.join([ variable, metric, 'cru_ts'+str(cru_ts_version), 'anom', month, year ])+'.tif' )
								for month, year in month_year ]

		# make a list of args to pass to the interpolation function
		args_list = [ {'anomalies':anom_df, 'meshgrid_tuple':(xi, yi), 'lons_pcll':lons_pcll, \
					'template_raster_fn':template_raster_fn, 'src_transform':src_transform, \
					'src_crs':src_crs, 'src_nodata':src_nodata, 'output_filename':fn } \
						for anom_df, fn in zip( anom_df_list, output_filenames ) ]
		
		anomalies = mp_map( lambda x: interpolate_anomalies( *x ), args_list, nproc=ncores )

		# read in the pre-processed 12-month climatology
		l = sorted( glob.glob( os.path.join( cl20_path, '*.tif' ) ) ) # this could catch you.
		clim_dict = { month:rasterio.open( fn ).read( 1 ) for month, fn in zip( months, l ) }

		# group the data by months
		out = pd.Series( out )
		out_months = out.apply( fn_month_grouper )
		months_grouped = out.groupby( out_months )

		# unpack groups for parallelization
		mg = [(i,j) for i,j in months_grouped ]
		# make an args tuple to pass to the function
		args_list = [ ( i[1], clim_dict[i[0]], downscaled_path, absolute ) for i in mg ]

		# downscale / write to disk
		out = mp_map( lambda args: downscale_cru_historical( *args ), args_list, nproc=ncores )
		return 'downscaling complete. files output at: %s' % base_path
	def downscale_cru_ts( cru_ts, clim_path, template_raster_fn, base_path, climatology_begin='1961', climatology_end='1990', ncores=2, absolute=True, metric='metric', variable=None, post_downscale_function=None ):
		'''
		downscale the cru TS data to a climatology template raster
		'''
		return main( **locals() )


if __name__ == '__main__':

	# this is a function that will perform a post-downscaling but pre-write modification of the interpolated
	# downscaled values, incase we got some values out of range while performing the interpolation with a
	# spline
	def clamp_vals( x ):
		''' clamp the values following the relative humidity downscaling '''
		x[ (x > 100) & (x < 500) ] = 95
		return x

	# run the downscaling of CRU TS3.*
	_ = DownscaleCRU.downscale_cru_ts( cru_ts, clim_path, template_raster_fn, base_path, climatology_begin='1961', climatology_end='1990', ncores=2, absolute=True, metric='metric', variable=None, post_downscale_function=None )



# # # # # HOW TO RUN THE APPLICATION # # # # # # # 
# # input args -- argparse it
# import os
# os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# ncores = '10'
# base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_ts31'
# cru_ts31 = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.cld.dat.nc' # 'cru_ts_3_10.1901.2009.tmp.nc'
# cl20_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
# template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
# anomalies_calc_type = 'relative' # 'absolute'
# downscaling_operation = 'mult' # 'add'

# climatology_begin = '1961'
# climatology_end = '1990'
# year_begin = '1901'
# year_end = '2009'
# variable = 'cld' # 'tas'
# metric = 'pct' # 'C'

# args_tuples = [ ('hi', cru_ts31), ('ci', cl20_path), ('tr', template_raster_fn), 
# 				('base', base_path), ('bt', year_begin), ('et', year_end),
# 				('cbt', climatology_begin), ('cet', climatology_end),
# 				('nc', ncores), ('at', anomalies_calc_type), ('m', metric), 
# 				('dso', downscaling_operation), ('v', variable) ]

# args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
# os.system( 'ipython2.7 -- tas_cld_cru_ts31_to_cl20_downscaling.py ' + args )

# # # # #TAS# # # # # # # 
# # input args -- argparse it
# import os
# os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# ncores = '5'
# base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_ts31'
# cru_ts31 = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.tmp.nc'
# cl20_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/tas/akcan'
# template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
# anomalies_calc_type = 'absolute'
# downscaling_operation = 'add'

# climatology_begin = '1961'
# climatology_end = '1990'
# year_begin = '1901'
# year_end = '2009'
# variable = 'tas'
# metric = 'C'

# args_tuples = [ ('hi', cru_ts31), ('ci', cl20_path), ('tr', template_raster_fn), 
# 				('base', base_path), ('bt', year_begin), ('et', year_end),
# 				('cbt', climatology_begin), ('cet', climatology_end),
# 				('nc', ncores), ('at', anomalies_calc_type), ('m', metric), 
# 				('dso', downscaling_operation), ('v', variable) ]

# args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
# os.system( 'ipython2.7 -- tas_cld_cru_ts31_to_cl20_downscaling.py ' + args )


# # # CRU TS 3.23 -- update:
# import os
# os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# ncores = '10'
# base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323'
# cru_ts31 = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323/cru_ts3.23.1901.2014.cld.dat.nc' # 'cru_ts_3_10.1901.2009.tmp.nc'
# cl20_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
# template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
# anomalies_calc_type = 'relative' # 'absolute'
# downscaling_operation = 'mult' # 'add'

# climatology_begin = '1961'
# climatology_end = '1990'
# year_begin = '1901'
# year_end = '2014'
# variable = 'cld' # 'tas'
# metric = 'pct' # 'C'

# args_tuples = [ ('hi', cru_ts31), ('ci', cl20_path), ('tr', template_raster_fn), 
# 				('base', base_path), ('bt', year_begin), ('et', year_end),
# 				('cbt', climatology_begin), ('cet', climatology_end),
# 				('nc', ncores), ('at', anomalies_calc_type), ('m', metric), 
# 				('dso', downscaling_operation), ('v', variable) ]

# args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
# os.system( 'ipython2.7 -- tas_cld_cru_ts31_to_cl20_downscaling.py ' + args )

# # TAS 
# import os
# os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# ncores = '10'
# base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323'
# cru_ts31 = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323/cru_ts3.23.1901.2014.tmp.dat.nc'
# cl20_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
# template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
# anomalies_calc_type = 'absolute'
# downscaling_operation = 'add'

# climatology_begin = '1961'
# climatology_end = '1990'
# year_begin = '1901'
# year_end = '2014'
# variable = tas'
# metric = 'C'

# args_tuples = [ ('hi', cru_ts31), ('ci', cl20_path), ('tr', template_raster_fn), 
# 				('base', base_path), ('bt', year_begin), ('et', year_end),
# 				('cbt', climatology_begin), ('cet', climatology_end),
# 				('nc', ncores), ('at', anomalies_calc_type), ('m', metric), 
# 				('dso', downscaling_operation), ('v', variable) ]

# args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
# os.system( 'ipython2.7 -- tas_cld_cru_ts31_to_cl20_downscaling.py ' + args )


