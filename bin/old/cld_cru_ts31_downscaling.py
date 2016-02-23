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
def shiftgrid(lon0,datain,lonsin,start=True,cyclic=360.0):
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
	template_meta = template_raster.meta

	if 'transform' in template_meta.keys():
		template_meta.pop( 'transform' )
		template_meta.update( crs={'init':'epsg:3338'} )

	reproject( dat, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
				dst_transform=template_meta['affine'], dst_crs=template_meta['crs'],\
				dst_nodata=None, resampling=RESAMPLING.nearest, num_threads=1, SOURCE_EXTRA=1000 )	
	return write_gtiff( output_arr, template_meta, output_filename, compress=True )

if __name__ == '__main__':
	import rasterio, xray, os, glob, affine
	from rasterio.warp import reproject, RESAMPLING
	import geopandas as gpd
	import pandas as pd
	import numpy as np
	from collections import OrderedDict
	from shapely.geometry import Point
	from pathos import multiprocessing as mp

	ncores = 15
	
	# filenames
	base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data'
	cld_ts31 = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.cld.dat.nc'
	template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	output_path = os.path.join( base_path, 'OCTOBER' )

	# this is the set of modified GTiffs produced in the conversion procedure with the ts2.0 data
	cld_ts20 = '' # read in the already pre-produced files.  They should be in 10'...  or maybe I need to change that.
	climatology_begin = '1961'
	climatology_end = '1990'
	year_begin = 1901
	year_end = 2009

	# open with xray
	cld_ts31 = xray.open_dataset( cld_ts31 )
	
	# open template raster
	template_raster = rasterio.open( template_raster_fn )
	template_meta = template_raster.meta
	template_meta.update( crs={'init':'epsg:3338'} )

	# calculate the anomalies
	clim_ds = cld_ts31.loc[ {'time':slice(climatology_begin,climatology_end)} ]
	climatology = clim_ds.cld.groupby( 'time.month' ).mean( 'time' )
	anomalies = cld_ts31.cld.groupby( 'time.month' ) / climatology

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
	years = np.arange( year_begin, year_end+1, 1 ).astype( str ).tolist()
	months = [ i if len(i)==2 else '0'+i for i in np.arange( 1, 12+1, 1 ).astype( str ).tolist() ]
	month_year = [ (month, year) for year in years for month in months ]
	output_filenames = [ os.path.join( output_path, '_'.join([ 'cld_pct_cru_ts31',month,year])+'.tif' ) for month, year in month_year ]

	# build a list of keyword args to pass to the pool of workers.
	args_list = [ {'df':df, 'meshgrid_tuple':(xi, yi), 'lons_pcll':lons_pcll, \
					'template_raster_fn':template_raster_fn, 'src_transform':src_transform, \
					'src_crs':src_crs, 'src_nodata':src_nodata, 'output_filename':fn } for df, fn in zip( df_list, output_filenames ) ]

	# interpolate / reproject / resample 
	pool = mp.Pool( processes=ncores )
	out = pool.map( lambda args: run( **args ), args_list )
	pool.close()

	# To Complete the CRU TS3.1 Downscaling we need the following: 
	# [1] DOWNSCALE WITH THE CRU CL2.0 Calculated Cloud Climatology from Sunshine Percent
	# [2] Mask the data 
	# [3] give proper naming convention
	# [4] output to GTiff
