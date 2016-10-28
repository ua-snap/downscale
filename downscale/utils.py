# -*- coding: utf8 -*-
# # # # 
# the beginnings of a class that will store utiity functions for downscaling
# these will be called by the individual downscaling scripts or read into the object
# as an object comprehension
# # # # 
import numpy as np
import rasterio

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
		See Rasterio documentation for more details. 

	RETURNS:
	--------
	string path to the new output_filename created

	'''
	import os, rasterio

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
	i0_shift = len(lonsin)-i0 # [ML] THIS COULD BE THE LINE GETTING US!!!
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
# def xyz_to_grid( x, y, z, grid, method='cubic', output_dtype=np.float32, *args, **kwargs ):
# 	'''
# 	interpolate points to a grid. simple wrapper around
# 	scipy.interpolate.griddata. Points and grid must be
# 	in the same coordinate system
# 	x = 1-D np.array of x coordinates / x,y,z must be same length
# 	y = 1-D np.array of y coordinates / x,y,z must be same length
# 	z = 1-D np.array of z coordinates / x,y,z must be same length
# 	grid = tuple of meshgrid as made using numpy.meshgrid()
# 			order (xi, yi)
# 	method = one of 'cubic', 'near', 'linear'
# 	'''
# 	from scipy.interpolate import griddata
# 	zi = griddata( (x, y), z, grid, method=method )
# 	zi = np.flipud( zi ).astype( output_dtype )
# 	return zi
def xyz_to_grid( x, y, z, grid, method='linear', output_dtype=np.float32, *args, **kwargs ):
	'''
	interpolate points to a grid. simple wrapper around
	matplotlib.mlab.griddata. Points and grid must be
	in the same coordinate system
	
	x = 1-D np.array of x coordinates / x,y,z must be same length
	y = 1-D np.array of y coordinates / x,y,z must be same length
	z = 1-D np.array of z coordinates / x,y,z must be same length
	grid = tuple of meshgrid as made using numpy.meshgrid()
			order (xi, yi)
	method = 'linear' -- hardwired currently and this is acceptable for
			a simple fill.
	'''
	from matplotlib.mlab import griddata
	xi, yi = grid
	zi = griddata( x, y, z, xi, yi, interp=method )
	return zi.astype( output_dtype )

def interp_ds( anom, base, src_crs, src_nodata, dst_nodata, src_transform, resample_type='bilinear',*args, **kwargs ):
	'''	
	anom = [numpy.ndarray] 2-d array representing a single monthly timestep of the data to be downscaled. 
							Must also be representative of anomalies.
	base = [str] filename of the corresponding baseline monthly file to use as template and downscale 
							baseline for combining with anomalies.
	src_transform = [affine.affine] 6 element affine transform of the input anomalies. [should be greenwich-centered]
	resample_type = [str] one of ['bilinear', 'count', 'nearest', 'mode', 'cubic', 'index', 'average', 'lanczos', 'cubic_spline']
	'''	
	import rasterio
	from rasterio.warp import reproject, RESAMPLING

	resampling = {'average':RESAMPLING.average,
				'cubic':RESAMPLING.cubic,
				'lanczos':RESAMPLING.lanczos,
				'bilinear':RESAMPLING.bilinear,
				'cubic_spline':RESAMPLING.cubic_spline,
				'mode':RESAMPLING.mode,
				'count':RESAMPLING.count,
				'index':RESAMPLING.index,
				'nearest':RESAMPLING.nearest }
	
	base = rasterio.open( base )
	baseline_arr = base.read( 1 )
	baseline_meta = base.meta
	baseline_meta.update( compress='lzw' )
	output_arr = np.empty_like( baseline_arr )
	
	reproject( anom, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
			dst_transform=baseline_meta['affine'], dst_crs=baseline_meta['crs'],\
			dst_nodata=dst_nodata, resampling=resampling[ resample_type ], SOURCE_EXTRA=1000 )
	return output_arr

def add( base, anom ):
	''' add anomalies to baseline '''
	return base + anom

def mult( base, anom ):
	''' multiply anomalies to baseline '''
	return base * anom

def _run_ds( d, f, operation_switch, anom=False, mask_value=0 ):
	'''
	run the meat of downscaling with this runner function for parallel processing

	ARGUMENTS:
	----------
	d = [dict] kwargs dict of args to pass to interpolation function
	f = [ ]

	RETURNS:
	--------

	'''
	import copy
	import rasterio
	
	post_downscale_function = d[ 'post_downscale_function' ]
	interped = f( **d )
	base = rasterio.open( d[ 'base' ] )
	base_arr = base.read( 1 )
	mask = base.read_masks( 1 )

	# set up output file metadata.
	meta = base.meta
	meta.update( compress='lzw' )
	if 'transform' in meta.keys():
		meta.pop( 'transform' )

	# write out the anomalies
	if anom == True:
		anom_filename = copy.copy( d[ 'output_filename' ] )
		dirname, basename = os.path.split( anom_filename )
		dirname = os.path.join( dirname, 'anom' )
		basename = basename.replace( '.tif', '_anom.tif' )
		try:
			if not os.path.exists( dirname ):
				os.makedirs( dirname )
		except:
			pass
		anom_filename = os.path.join( dirname, basename )
		with rasterio.open( anom_filename, 'w', **meta ) as anom:
			anom.write( interped, 1 )
	
	# make sure the output dir exists and if not, create it
	dirname = os.path.dirname( d[ 'output_filename' ] )
	if not os.path.exists( dirname ):
		os.makedirs( dirname )

	# operation switch
	output_arr = operation_switch[ d[ 'downscaling_operation' ] ]( base_arr, interped )
	
	# post downscale it if func given
	if post_downscale_function != None:
		output_arr = post_downscale_function( output_arr )
		# drop the mask if there is one
		if hasattr( output_arr, 'mask'):
			output_arr = output_arr.data

	# make sure data is masked
	output_arr[ mask == mask_value ] = meta[ 'nodata' ]

	# write it to disk.
	with rasterio.open( d[ 'output_filename' ], 'w', **meta ) as out:
		out.write( output_arr, 1 )
	return d['output_filename']


# def downscale( anom_arr, baseline_arr, output_filename,	downscaling_operation, \
# 	meta, post_downscale_function, mask=None, mask_value=0, *args, **kwargs ):
# 	'''
# 	downscale an anomaly array with a baseline array from the same period.

# 	Arguments:
# 	----------
# 	anom_arr = [ np.ndarray ] 2-D NumPy array representing a raster domain. 
# 				anom/baseline arrays must be same shape.
# 	baseline_arr = [ np.ndarray ] 2-D NumPy array representing a raster domain. 
# 				anom/baseline arrays must be same shape.
# 	output_filename = [ str ] full path and output filename to be created
# 	downscaling_operation = [ ] 
# 	meta = [ dict ] rasterio-style dictionary of raster metadata attributes. This 
# 			must jive with the dimensions and the data type of the array generated 
# 			through downscaling anom_arr with baseline_arr.  
# 	post_downscale_function = [ function ] a function that takes a 2-D downscaled 
# 			array as input and returns an array of the same shape / datatype.  This
# 			is typically used as a post-mortem for clamping the values from an output
# 			downscaled array that may be slightly outside the range due to the 
# 			interpolation method. We currently use this to clamp the values of the hur
# 			to 0-100.

# 	Returns:
# 	--------
# 	output_filename of newly generated downscaled raster.

# 	'''
# 	import rasterio

# 	def add( base, anom ):
# 		return base + anom
# 	def mult( base, anom ):
# 		return base * anom
# 	def div( base, anom ):
# 		# this one may not be useful, but the placeholder is here
# 		# return base / anom
# 		return NotImplementedError

# 	try:
# 		operation_switch = { 'add':add, 'mult':mult, 'div':div }
# 	except:
# 		AttributeError( 'downscale: incorrect downscaling_operation str' )
	
# 	output_arr = operation_switch[ downscaling_operation ]( baseline_arr, anom_arr )
# 	output_arr[ np.isinf( output_arr ) ] = meta[ 'nodata' ]

# 	if isinstance(mask, np.ndarray):
# 		output_arr[ mask == 0 ] = mask_value

# 	if post_downscale_function != None:
# 		output_arr = post_downscale_function( output_arr )

# 	if 'transform' in meta.keys():
# 		# avoid the gdal geotransform deprecation warning
# 		meta.pop( 'transform' )

# 	with rasterio.open( output_filename, 'w', **meta ) as out:
# 		out.write( output_arr, 1 )
# 	return output_filename


# OLDER NON-SNAPPING version:
# def transform_from_latlon( lat, lon ):
# 	''' simple way to make an affine transform from lats and lons coords '''
# 	from affine import Affine
# 	lat = np.asarray( lat )
# 	lon = np.asarray( lon )
# 	trans = Affine.translation(lon[0], lat[0])
# 	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
# 	return trans * scale

def transform_from_latlon( lat, lon ):
	''' simple way to make an affine transform from lats and lons coords '''
	from affine import Affine
	lat = np.asarray( lat )
	lon = np.asarray( lon )
	if (np.max( lat ) - 90) < np.abs( np.mean( np.diff( lat ) ) ):
		lat_max = 90.0
	else:
		lat_max = np.max( lat )

	# set the lonmax to the corner. --> this can get you into trouble with non-global data
	# but I am unsure how to make it more dynamic at the moment. [ML]
	lon_arr = np.array([-180.0, 0.0 ])
	idx = (np.abs(lon_arr - np.min( lon ) ) ).argmin()
	lon_max = lon_arr[ idx ]

	trans = Affine.translation(lon_max, lat_max)
	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
	return trans * scale

def rasterize( shapes, coords, latitude='lat', longitude='lon', fill=None, **kwargs ):
	'''
	Rasterize a list of (geometry, fill_value) tuples onto the given
	xarray coordinates. This only works for 1d latitude and longitude
	arrays.

	ARGUMENTS:
	----------
	shapes = [list] of tuples of (shapely.geom, fill_value)
	coords = [dict] of named 1d latitude and longitude arrays.
	latitude = [str] name of latitude key. default:'latitude'
	longitude = [str] name of longitude key. default:'longitude'
	fill = fill_value

	RETURNS:
	--------

	xarray.DataArray

	'''
	from rasterio import features
	import xarray as xr
	if fill == None:
		fill = np.nan
	transform = transform_from_latlon( coords[ latitude ], coords[ longitude ] )
	out_shape = ( len( coords[ latitude ] ), len( coords[ longitude ] ) )
	raster = features.rasterize( shapes, out_shape=out_shape,
								fill=fill, transform=transform,
								dtype=float, **kwargs )
	spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
	return xr.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))


# # EXAMPLE OF HOW TO RASTERIZE A SHAPE TO An xarray.Dataset
# # this shapefile is from natural earth data
# # http://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/
# states = geopandas.read_file('/Users/shoyer/Downloads/ne_10m_admin_1_states_provinces_lakes')
# us_states = states.query("admin == 'United States of America'").reset_index(drop=True)
# state_ids = {k: i for i, k in enumerate(us_states.woe_name)}
# shapes = [(shape, n) for n, shape in enumerate(us_states.geometry)]

# ds = xray.Dataset(coords={'longitude': np.linspace(-125, -65, num=5000),
#                           'latitude': np.linspace(50, 25, num=3000)})
# ds['states'] = rasterize(shapes, ds.coords)

# example of applying a mask
# ds.states.where(ds.states == state_ids['California'])

# # # # # # # # # END! NEW FILL Dataset
