# -*- coding: utf8 -*-
# # # # 
# the beginnings of a class that will store utiity functions for downscaling
# these will be called by the individual downscaling scripts or read into the object
# as an object comprehension
# # # # 

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
def xyz_to_grid( x, y, z, grid, method='cubic', output_dtype=np.float32, *args, **kwargs  ):
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
def downscale( anom_arr, baseline_arr, output_filename, \
	downscaling_operation, meta, post_downscale_function, mask=None, mask_value=0, *args, **kwargs ):
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
	import rasterio

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
	
	output_arr = operation_switch[ downscaling_operation ]( baseline_arr, anom_arr )
	output_arr[ np.isinf( output_arr ) ] = meta[ 'nodata' ]

	if isinstance(mask, np.ndarray):
		output_arr[ mask == 1 ] = mask_value

	if post_downscale_function != None:
		output_arr = post_downscale_function( output_arr )

	if 'transform' in meta.keys():
		# avoid the gdal geotransform deprecation warning
		meta.pop( 'transform' )

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( output_arr, 1 )
	return output_filename
