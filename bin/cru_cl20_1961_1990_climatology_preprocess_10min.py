# THIS SCRIPT WILL convert the point csv-like file from the CRU CL2.0 and will interpolate it to 
# a worldwide grid at 10' resolution along with an akcan grid (matching SNAP baseline) at 2km resolution.
# 
# written by: Michael Lindgren (malindgren@alaska.edu)
# The MIT License (MIT)
# # # # # # 
# Copyright (c) 2015 Michael A. Lindgren
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# # # # # # 

import numpy as np # hack to solve a lib issue in the function args of xyztogrid
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
	if np.fabs( lonsin[-1]-lonsin[0]-cyclic ) > 1.e-4:
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
def cru_xyz_to_shp( in_xyz, lon_col, lat_col, crs, output_filename=None ):
	'''
	convert the cru cl2.0 1961-1990 Climatology data to a shapefile.
		*can handle the .dat format even if compressed with .gzip extension.
	PARAMETERS:
	-----------
	in_xyz = path to the .dat or .dat.gz downloaded cru cl2.0 file from UK Met Office site
	lon_col = string name of column storing longitudes
	lat_col = string name of column storing latitudes
	crs = proj4string or epsg code
	output_filename = string path to the output filename to be created
						*if None it will return only the new GeoDataFrame
	RETURNS
	-------
	output_filename as string or GeoDataFrame depending on the value of 
		output_filename
	'''
	colnames = ['lat', 'lon', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
	from shapely.geometry import Point
	import pandas as pd
	import geopandas as gpd
	import os

	if os.path.splitext( in_xyz )[1] == '.gz':
		cru_df = pd.read_csv( in_xyz, delim_whitespace=True, compression='gzip', header=None, names=colnames )
	else:
		cru_df = pd.read_csv( in_xyz, delim_whitespace=True, header=None, names=colnames )

	# create a column named geometry with shapely geometry objects for each row
	def f( x ):
		''' return a Point shapely object for each x,y pair'''
		return Point( x[1], x[0] ) # Point( x.lon, x.lat ) # HACK!
	
	cru_df[ 'geometry' ] = cru_df.apply( f, axis=1 )
	cru_df = gpd.GeoDataFrame( cru_df ) # convert to GeoDataFrame
	if output_filename:
		cru_df.to_file( output_filename, 'ESRI Shapefile' )
		return output_filename
	else:
		return cru_df
def bounds_to_extent( bounds ):
	'''
	take input rasterio bounds object and return an extent
	'''
	l,b,r,t = bounds
	return [ (l,b), (r,b), (r,t), (l,t), (l,b) ]
def extent_to_shapefile( extent, output_shapefile, proj4string ):
	''' convert an extent to a shapefile using its proj4string '''
	import geopandas as gpd
	from shapely.geometry import Polygon
	gpd.GeoDataFrame( {'extent_id':1, 'geometry':Polygon( extent )}, index=[1], crs=proj4string ).to_file( output_shapefile, 'ESRI Shapefile' )
	return output_shapefile
def pad_bounds( rst, npixels, crs, output_shapefile ):
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
	
	new_ext = bounds_to_extent( new_bounds )
	return extent_to_shapefile( new_ext, output_shapefile, crs )
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
	method = one of 'cubic', 'near', 'linear'
	'''
	import numpy as np
	from scipy.interpolate import griddata
	zi = griddata( (x, y), z, grid, method=method )
	return zi
def interpolate_to_grid( x, y, z, grid, output_filename, meta, mask=None, method='cubic', output_dtype=np.float32 ):
	'''
	x: list of unique longitudes
	y: list of unique latitudes
	z: list of values for the grid
	grid: lons and lats mesh gridded
	output_filename: string path to the new output raster to create
	meta: rasterio-style metadata dictionary describing a raster
	mask: 2d array in the same shape as the grid, to be used to mask the output array before writing
			to disk.
	method: one of ['cubic', 'nearest', 'linear'] to use as the interpolation methodology
	output_dtype: output datatype (numpy)
	'''
	import numpy as np
	import rasterio, os
	# interpolate to the grid
	interp_arr = xyz_to_grid( x, y, z, grid, method=method, output_dtype=output_dtype )
	# write it to disk
	meta.update( compress='lzw' )
	with rasterio.open( output_filename, 'w', **meta ) as out:
		if mask:
			mask = template_rst.read_masks( 1 )
			interp_arr[ mask == 0 ] = meta[ 'nodata' ]
		out.write_band( 1, interp_arr.astype( output_dtype ) )
	return output_filename
def generate_10min_grid():
	from affine import Affine as A
	resolution = 0.166667
	l, b, r, t = [ -180.0, -90, 180.0, 90 ] # the global extent

	cols = float( round( ( r-l )/resolution ) )
	rows = float( round( ( t-b )/resolution ) )
	affine_transform = A( (360/cols), 0.0, l, 0.0, -(180/rows), t )
	crs = {'init':'epsg:4326'}

	meta = {}
	meta.update( compress='lzw', affine=affine_transform, width=cols, height=rows, crs=crs, nodata=None, dtype=np.float32, count=1, driver='GTiff' )
	# this is the new grid that we will interpolate the points to
	grid = np.zeros( shape=(rows,cols), dtype=np.float32 )
	return meta, grid
def coordinates( fn=None, meta=None, numpy_array=None, input_crs=None, to_latlong=False ):
	'''
	take a raster file as input and return the centroid coords for each 
	of the grid cells as a pair of numpy 2d arrays (longitude, latitude)
	'''
	import rasterio
	import numpy as np
	from affine import Affine
	from pyproj import Proj, transform

	if fn:
		# Read raster
		with rasterio.open( fn ) as r:
			T0 = r.affine  # upper-left pixel corner affine transform
			p1 = Proj( r.crs )
			A = r.read_band( 1 )  # pixel values

	elif (meta is not None) & (numpy_array is not None):
		A = numpy_array
		if input_crs != None:
			p1 = Proj( input_crs )
			T0 = meta[ 'affine' ]
		else:
			p1 = None
			T0 = meta[ 'affine' ]
	else:
		BaseException( 'check inputs' )

	# All rows and columns
	cols, rows = np.meshgrid(np.arange(A.shape[1]), np.arange(A.shape[0]))
	# Get affine transform for pixel centres
	T1 = T0 * Affine.translation( 0.5, 0.5 )
	# Function to convert pixel row/column index (from 0) to easting/northing at centre
	rc2en = lambda r, c: ( c, r ) * T1
	# All eastings and northings (there is probably a faster way to do this)
	eastings, northings = np.vectorize(rc2en, otypes=[np.float, np.float])(rows, cols)

	if to_latlong == False:
		return eastings, northings
	elif (to_latlong == True) & (input_crs != None):
		# Project all longitudes, latitudes
		longs, lats = transform(p1, p1.to_latlong(), eastings, northings)
		return longs, lats
	else:
		BaseException( 'cant reproject to latlong without an input_crs' )
def regrid( src, dst, src_affine, src_crs, dst_affine, dst_crs, method='cubic_spline' ):
	from rasterio.warp import RESAMPLING, reproject
	resampling_type = {'nearest':RESAMPLING.nearest,
					'bilinear':RESAMPLING.bilinear,
					'cubic':RESAMPLING.cubic,
					'cubic_spline':RESAMPLING.cubic_spline,
					'lanczos':RESAMPLING.lanczos,
					'average':RESAMPLING.average,
					'mode':RESAMPLING.mode}

	out = np.zeros( dst.shape )
	reproject( src,
				out,
				src_transform=src_affine,
				src_crs=src_crs,
				dst_transform=dst_affine,
				dst_crs=dst_crs,
				resampling=resampling_type[ method ], num_threads=2, SOURCE_EXTRA=5000 )
	return out
def run( args ):
	# unpack some args being passed into the run wrapper
	x = args['x']
	y = args['y']
	z = args['z']
	meshgrid_10min = args['meshgrid_10min']
	output_filename_10min = args['output_filename_10min']
	meta_10min = args['meta_10min']
	meta_10min.update( compress='lzw' )

	# interpolate to a global (unmasked) 10 min output
	new_grid = interpolate_to_grid( x, y, z, meshgrid_10min, output_filename_10min, meta_10min )
	new_grid = rasterio.open( output_filename_10min ).read( 1 ) # read it back in << this is a hack!

	# unpack some args being passed in
	template_raster = args['template_raster']
	src_affine = args['meta_10min']['affine']
	dst_affine = template_raster.affine
	method = 'cubic_spline'
	output_filename = args['output_filename']
	src_crs = {'init':'epsg:4326'}
	dst_crs = {'init':'epsg:3338'}
	meta_akcan = template_raster.meta
	meta_akcan.update( compress='lzw', crs=dst_crs  )
	mask = args['mask']

	# regrid to 3338 and AKCAN extent -- currently hardwired
	akcan = regrid( new_grid, template_raster.read( 1 ), src_affine, src_crs, dst_affine, dst_crs, method=method )
	with rasterio.open( output_filename, 'w', **meta_akcan ) as out:
		if mask is not None:
			akcan[ mask == 0 ] = meta_akcan[ 'nodata' ]
		out.write( akcan.astype( np.float32 ), 1 ) # watch this hardwired type!
	return output_filename

if __name__ == '__main__':
	import os, rasterio, glob, fiona
	import scipy as sp
	import numpy as np
	import pandas as pd
	import geopandas as gpd
	from rasterio import Affine as A
	from pathos import multiprocessing as mp
	from shapely.geometry import Point
	import argparse

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='preprocess CRU CL2.0 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-p", "--base_path", action='store', dest='base_path', type=str, help="path to parent directory with a subdirector(ies)y storing the data" )
	parser.add_argument( "-cru", "--cru_filename", action='store', dest='cru_filename', type=str, help="string path to the .tar.gz file location, downloaded from the CRU site" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="string abbreviated name of the variable being processed." )
	parser.add_argument( "-tr", "--template_raster_fn", action='store', dest='template_raster_fn', type=str, help="string path to a template raster dataset to match the CRU CL2.0 to." )

	# parse and unpack the args
	args = parser.parse_args()
	base_path = args.base_path
	cru_filename = args.cru_filename
	variable = args.variable
	template_raster_fn = args.template_raster_fn
	
	# read in the template output raster
	template_raster = rasterio.open( template_raster_fn )

	# make a 10min global grid to store the newly interpolated (to overcome coast issues) CL2.0 files
	meta_10min, grid_10min = generate_10min_grid() # NEW! Greenwich centered...
	# update the metadata to match the 0-360 longitude ordering
	meta_10min.update( affine=A(0.16666666666666666, 0.0, 0, 0.0, -0.16666666666666666, 90.0)) # pacific centered...

	# build an output path to store the data generated with this script
	cru_path = os.path.join( base_path, 'cru_ts20', variable )

	if not os.path.exists( cru_path ):
		os.makedirs( cru_path )

	# convert to point geopandas geodataframe
	# cru_shp_fn = os.path.join( cru_path, 'cru_'+variable+'_ts20_1961_1990_climatology.shp' )
	cru_gdf = cru_xyz_to_shp( cru_filename, 'lon', 'lat', {'init':'epsg:4326'}, output_filename=None )
	cru_gdf.loc[ cru_gdf.lon < 0, 'lon' ] = cru_gdf.loc[ cru_gdf.lon < 0, 'lon' ] + 360.0 # this updates the lon to PCLL
	
	# # # # [ TEST ] WRITE OUT THE SHAPEFILE # # # # <- this works
	# cru_gdf[ 'geometry' ] = cru_gdf[ ['lon', 'lat'] ].apply( lambda x: Point(*x),  axis=1 )
	# cru_gdf.to_file( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_v4/test_pcll_shape.shp', 'ESRI Shapefile' )
	# # # # [ END TEST ] WRITE OUT THE SHAPEFILE # # # #

	# get the lons and lats in a grid and then rotate them to a 0-360 -- Pacific Centered Latlong
	lons, lats = coordinates( meta=meta_10min, numpy_array=grid_10min )

	# build some args
	months = ['01','02','03','04','05','06','07','08','09','10','11','12']
	output_filenames_10min = [ os.path.join( cru_path, variable+'_cru_cl20_10min_'+month+'_1961_1990.tif' ) for month in months ]
	output_filenames = [ os.path.join( cru_path, variable+'_cru_cl20_akcan_'+month+'_1961_1990.tif' ) for month in months ]
	akcan_mask = template_raster.read_masks( 1 )

	# run in parallel
	args_list = [ { 'x':cru_gdf['lon'], 'y':cru_gdf['lat'], 'z':np.array(cru_gdf[ month ]), \
					'meshgrid_10min':(lons,lats), 'output_filename_10min':out_fn_10min, 'output_filename':out_fn, \
					'meta_10min':meta_10min, 'meta_akcan':template_raster.meta, 'mask':akcan_mask, 'method':'cubic', 'template_raster':template_raster } \
						for month, out_fn_10min, out_fn in zip( months, output_filenames_10min, output_filenames ) ]

	pool = mp.ThreadPool( 2 )
	out = pool.map( run, args_list )
	pool.close()
	# out = map( run, args_list )

# # # EXAMPLE OF USE # # #
# import os
# os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# cru_folder = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20'
# var_fn_dict = { 'hur':os.path.join( cru_folder, 'grid_10min_reh.dat.gz'),'tas':os.path.join( cru_folder, 'grid_10min_tmp.dat.gz'), 'sunp':os.path.join( cru_folder, 'grid_10min_sunp.dat.gz' ) }
# base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_v4'
# template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'

# for variable, cru_filename in var_fn_dict.iteritems():
# 	print 'working on : %s' % variable
# 	os.system( 'ipython cru_cl20_1961_1990_climatology_preprocess_10min.py -- -p ' + base_path + ' -cru ' + cru_filename + ' -v ' + variable + ' -tr ' + template_raster_fn )
