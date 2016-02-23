# # #
# Current implementation of the cru ts31 (ts32) delta downscaling procedure
#
# Author: Michael Lindgren (malindgren@alaska.edu)
# # #
import numpy as np
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
		supported here.

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
	# new_ext = bounds_to_extent( new_bounds )
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
	import numpy as np
	from scipy.interpolate import griddata

	zi = griddata( (x, y), z, grid, method=method )
	zi = np.flipud( zi.astype( output_dtype ) )
	return zi
def run( df, meshgrid_tuple, lons_pcll, template_raster_fn, src_transform, src_crs, src_nodata, output_filename ):
	'''
	run the interpolation to a grid, and reprojection / resampling to the Alaska / Canada rasters
	extent, resolution, origin (template_raster).

	This function is intended to be used to run a pathos.multiprocessing Pool's map function
	across a list of pre-computed arguments.
			
	RETURNS:

	[str] path to the output filename generated

	'''
	template_raster = rasterio.open( template_raster_fn )
	interp_arr = xyz_to_grid( np.array(df['lon'].tolist()), \
					np.array(df['lat'].tolist()), \
					np.array(df['anom'].tolist()), grid=meshgrid_tuple, method='cubic' ) 

	src_nodata = -9999.0 # nodata
	interp_arr[ np.isnan( interp_arr ) ] = src_nodata
	dat, lons = shiftgrid( 180., interp_arr, lons_pcll, start=False )
	output_arr = np.empty_like( template_raster.read( 1 ) )
	# mask it with the internal mask in the template raster, where 0 is oob.
	output_arr = np.ma.masked_where( template_raster.read_masks( 1 ) == 0, output_arr )
	template_meta = template_raster.meta

	if 'transform' in template_meta.keys():
		template_meta.pop( 'transform' )

	reproject( dat, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
				dst_transform=template_meta['affine'], dst_crs=template_meta['crs'],\
				dst_nodata=None, resampling=RESAMPLING.nearest, num_threads=1, SOURCE_EXTRA=1000 )	
	return write_gtiff( output_arr, template_meta, output_filename, compress=True )
def fn_month_grouper( x ):
	'''
	take a filename and return the month element of the naming convention
	'''
	return os.path.splitext(os.path.basename(x))[0].split( '_' )[5]
def convert_to_vap( tas_arr, hur_arr ):
	''' create relative humidity from the CRU tas / vap '''
	esa_arr = 6.112 * np.exp( 17.62 * tas_arr/ (243.12 + tas_arr) )
	# esa_arr = 6.112 * np.exp( 22.46 * tas_arr / (272.62 + tas_arr) )
	return (hur_arr*esa_arr)/100
def downscale_cru_historical( file_list, cru_cl20_arr, output_path, downscaling_operation ):
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
	
	def f( anomaly_fn, baseline_arr, output_path ):
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
		meta.update( compress='lzw' )
		cru_ts31 = cru_ts31.read( 1 )
		operation_switch = { 'add':add, 'mult':mult, 'div':div }
		downscaled = operation_switch[ downscaling_operation ]( baseline_arr, cru_ts31 )
		# this is hardwired stuff for this fairly hardwired script.
		output_filename = os.path.basename( anomaly_fn ).replace( 'anom', 'downscaled' )
		output_filename = os.path.join( output_path, output_filename )
		output_arr = baseline_arr * cru_ts31 # multiply since it was relative anomalies
		if 'transform' in meta.keys():
			meta.pop( 'transform' )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( output_arr, 1 )
		return output_filename
		
	partial_f = partial( f, baseline_arr=cru_cl20_arr, output_path=output_path )
	cru_ts31 = file_list.apply( lambda fn: partial_f( anomaly_fn=fn ) )
	return output_path


if __name__ == '__main__':
	import rasterio, xray, os, glob, affine
	from rasterio.warp import reproject, RESAMPLING
	import geopandas as gpd
	import pandas as pd
	import numpy as np
	from collections import OrderedDict
	from shapely.geometry import Point
	from pathos import multiprocessing as mp
	import argparse

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='preprocess cmip5 input netcdf files to a common type and single files' )
	parser.add_argument( "-hhi", "--cru_ts31_vap", action='store', dest='cru_ts31_vap', type=str, help="path to historical CRU TS3.1 vap input NetCDF file" )
	parser.add_argument( "-thi", "--cru_ts31_tas", action='store', dest='cru_ts31_tas', type=str, help="path to historical CRU TS3.1 tas input NetCDF file" )
	parser.add_argument( "-cit", "--cl20_tas_path", action='store', dest='cl20_tas_path', type=str, help="path to historical TAS (in output extent/res/crs) CRU TS2.0 Climatology input directory in single-band GTiff Format" )
	parser.add_argument( "-cih", "--cl20_hur_path", action='store', dest='cl20_hur_path', type=str, help="path to historical HUR (in output extent/res/crs) CRU TS2.0 Climatology input directory in single-band GTiff Format" )
	parser.add_argument( "-tr", "--template_raster_fn", action='store', dest='template_raster_fn', type=str, help="path to ALFRESCO Formatted template raster to match outputs to." )
	parser.add_argument( "-base", "--base_path", action='store', dest='base_path', type=str, help="string path to the folder to put the output files into" )
	parser.add_argument( "-bt", "--year_begin", action='store', dest='year_begin', type=int, help="string in format YYYY of the beginning year in the series" )
	parser.add_argument( "-et", "--year_end", action='store', dest='year_end', type=int, help="string in format YYYY of the ending year in the series" )
	parser.add_argument( "-cbt", "--climatology_begin_time", nargs='?', const='196101', action='store', dest='climatology_begin', type=str, help="string in format YYYY of the beginning year of the climatology period" )
	parser.add_argument( "-cet", "--climatology_end_time", nargs='?', const='199012', action='store', dest='climatology_end', type=str, help="string in format YYYY of the ending year of the climatology period" )
	parser.add_argument( "-nc", "--ncores", nargs='?', const=2, action='store', dest='ncores', type=int, help="integer valueof number of cores to use. default:2" )
	parser.add_argument( "-at", "--anomalies_calc_type", nargs='?', const='absolute', action='store', dest='anomalies_calc_type', type=str, help="string of 'proportional' or 'absolute' to inform of anomalies calculation type to perform." )
	parser.add_argument( "-m", "--metric", nargs='?', const='metric', action='store', dest='metric', type=str, help="string of whatever the metric type is of the outputs to put in the filename." )
	parser.add_argument( "-dso", "--downscaling_operation", action='store', dest='downscaling_operation', type=str, help="string of 'add', 'mult', 'div', which refers to the type or downscaling operation to use." )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="string of the abbreviation used to identify the variable (i.e. cld)." )

	# parse args
	args = parser.parse_args()

	# unpack args
	ncores = args.ncores
	base_path = args.base_path
	cru_ts31_vap = args.cru_ts31_vap
	cru_ts31_tas = args.cru_ts31_tas
	cl20_tas_path = args.cl20_tas_path
	cl20_hur_path = args.cl20_hur_path
	template_raster_fn = args.template_raster_fn
	anomalies_calc_type = args.anomalies_calc_type
	downscaling_operation = args.downscaling_operation
	climatology_begin = args.climatology_begin
	climatology_end = args.climatology_end
	year_begin = args.year_begin
	year_end = args.year_end
	variable = args.variable
	metric = args.metric

	# make some output directories if they are not there already to dump 
	# our output files
	anomalies_path = os.path.join( base_path, variable, 'anom' )
	if not os.path.exists( anomalies_path ):
		os.makedirs( anomalies_path )

	downscaled_path = os.path.join( base_path, variable, 'downscaled' )
	if not os.path.exists( downscaled_path ):
		os.makedirs( downscaled_path )

	# open with xray
	cru_ts31_vap = xray.open_dataset( cru_ts31_vap )
	
	# open template raster
	template_raster = rasterio.open( template_raster_fn )
	template_meta = template_raster.meta
	template_meta.update( crs={'init':'epsg:3338'} )

	# make a mask with values of 0=nodata and 1=data
	template_raster_mask = template_raster.read_masks( 1 )
	template_raster_mask[ template_raster_mask == 255 ] = 1

	# calculate the anomalies
	# this is temporary name change for the tmp (tas) data naming diff.
	if variable == 'tas':
		variable = 'tmp'
	
	clim_ds = cru_ts31.loc[ {'time':slice(climatology_begin,climatology_end)} ]
	climatology = clim_ds[ variable ].groupby( 'time.month' ).mean( 'time' )

	if anomalies_calc_type == 'relative':
		anomalies = cru_ts31[ variable ].groupby( 'time.month' ) / climatology

	if anomalies_calc_type == 'absolute':
		anomalies = cru_ts31[ variable ].groupby( 'time.month' ) - climatology

	# reset the variable if tas
	if variable == 'tmp':
		variable = 'tas'

	# rotate the anomalies to pacific centered latlong -- this is already in the greenwich latlong
	dat_pcll, lons_pcll = shiftgrid( 0., anomalies, anomalies.lon.data )
	
	# # generate an expanded extent (from the template_raster) to interpolate across
	template_raster = rasterio.open( template_raster_fn )
	# output_resolution = (1000.0, 1000.0) # hardwired, but we are building this for IEM which requires 1km
	template_meta = template_raster.meta

	# # interpolate to a new grid
	# get longitudes and latitudes using meshgrid
	lo, la = [ i.ravel() for i in np.meshgrid( lons_pcll, anomalies.lat ) ] # mesh the lons/lats
	
	# convert into GeoDataFrame and drop all the NaNs
	df_list = [ pd.DataFrame({ 'anom':i.ravel(), 'lat':la, 'lon':lo }).dropna( axis=0, how='any' ) for i in dat_pcll ]
	xi, yi = np.meshgrid( lons_pcll, anomalies.lat.data )
	# meshgrid_tuple = np.meshgrid( lons_pcll, anomalies.lat.data )

	# argument setup
	src_transform = affine.Affine( 0.5, 0.0, -180.0, 0.0, -0.5, 90.0 )
	src_crs = {'init':'epsg:4326'}
	src_nodata = -9999.0
		
	# output_filenames setup
	years = np.arange( int(year_begin), int(year_end)+1, 1 ).astype( str ).tolist()
	months = [ i if len(i)==2 else '0'+i for i in np.arange( 1, 12+1, 1 ).astype( str ).tolist() ]
	month_year = [ (month, year) for year in years for month in months ]
	output_filenames = [ os.path.join( anomalies_path, '_'.join([ variable,metric,'cru_ts31_anom',month,year])+'.tif' ) 
							for month, year in month_year ]

	# build a list of keyword args to pass to the pool of workers.
	args_list = [ {'df':df, 'meshgrid_tuple':(xi, yi), 'lons_pcll':lons_pcll, \
					'template_raster_fn':template_raster_fn, 'src_transform':src_transform, \
					'src_crs':src_crs, 'src_nodata':src_nodata, 'output_filename':fn } \
						for df, fn in zip( df_list, output_filenames ) ]

	# interpolate / reproject / resample the anomalies to match template_raster_fn
	pool = mp.Pool( processes=ncores )
	out = pool.map( lambda args: run( **args ), args_list )
	pool.close()

	# To Complete the CRU TS3.1 Downscaling we need the following:
	# read in the pre-processed CL2.0 Cloud Climatology
	# open hur / tas CRU TS2.0 Climatologies and convert the hur to vap
	tas_clim = glob.glob( os.path.join( cl20_tas_path, '*.tif' ) )
	hur_clim = glob.glob( os.path.join( cl20_hur_path, '*.tif'  ) )

	cl20 = { month:convert_to_vap( rasterio.open(tas_fn).read(1), rasterio.open(hur_fn).read(1) ) \
				for month, tas_fn, hur_fn in zip( months, tas_clim, hur_clim) }

	l = sorted( glob.glob( os.path.join( cl20_path, '*.tif' ) ) ) # this could catch you.
	cl20_dict = { month:rasterio.open( fn ).read( 1 ) for month, fn in zip( months, l ) }

	# group the data by months
	out = pd.Series( out )
	out_months = out.apply( fn_month_grouper )
	months_grouped = out.groupby( out_months )

	# unpack groups for parallelization and make a list of tuples of arguments to pass to the downscale function
	mg = [(i,j) for i,j in months_grouped ]
	args_list = [ ( i[1], cl20_dict[i[0]], downscaled_path, downscaling_operation ) for i in mg ]

	# downscale / write to disk
	pool = mp.Pool( processes=ncores )
	out = pool.map( lambda args: downscale_cru_historical( *args ), args_list )
	pool.close()


# # # # # HOW TO RUN THE APPLICATION # # # # # # # 
# # input args -- argparse it
# import os
# os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# ncores = '10'
# base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_ts31'
# cru_ts31_vap = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.vap.nc'
# cru_ts31_tas = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.tmp.nc'
# cl20_tas_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_v2/cru_ts20/tas/akcan'
# cl20_hur_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_v2/cru_ts20/hur/akcan'
# template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
# anomalies_calc_type = 'absolute' # 'relative'
# downscaling_operation = 'add' # 'mult', 'div'

# climatology_begin = '1961'
# climatology_end = '1990'
# year_begin = '1901'
# year_end = '2009'
# variable = 'vap'
# metric = 'hPa'

# args_tuples = [ ('hhi', cru_ts31_vap), ('thi', cru_ts31_tas), ('ci', cl20_path), ('tr', template_raster_fn), 
# 				('base', base_path), ('bt', year_begin), ('et', year_end), 
# 				('cbt', climatology_begin), ('cet', climatology_end), 
# 				('nc', ncores), ('at', anomalies_calc_type), ('m', metric), 
# 				('dso', downscaling_operation), ('v', variable) ]

# args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
# os.system( 'python cru_ts31_to_cl20_downscaling.py ' + args )


