# preprocess prism using points to a known grid...
# UNUSED METHOD BUILT TO TEST HOW OLD DATA WAS MADE.  THIS METHOD LEAVES ARTIFACTS.
import numpy as np
def coordinates( fn=None, meta=None, numpy_array=None, input_crs=None, to_latlong=False ):
	'''
	take a raster file as input and return the centroid coords for each 
	of the grid cells as a pair of numpy 2d arrays (longitude, latitude)
	User must give either:
		fn = path to the rasterio readable raster
	OR
		meta & numpy ndarray (usually obtained by rasterio.open(fn).read( 1 )) 
		where:
		meta = a rasterio style metadata dictionary ( rasterio.open(fn).meta )
		numpy_array = 2d numpy array representing a raster described by the meta
	input_crs = rasterio style proj4 dict, example: { 'init':'epsg:3338' }
	to_latlong = boolean.  If True all coordinates will be returned as EPSG:4326
						 If False all coordinates will be returned in input_crs
	returns:
		meshgrid of longitudes and latitudes
	''' 
	
	import rasterio
	import numpy as np
	from affine import Affine
	from pyproj import Proj, transform

	if fn:
		# Read raster
		with rasterio.open( fn ) as r:
			T0 = r.transform  # upper-left pixel corner affine transform
			p1 = r.crs
			A = r.read( 1 )  # pixel values

	elif (meta is not None) & (numpy_array is not None):
		A = numpy_array
		if input_crs != None:
			p1 = input_crs
			T0 = meta[ 'transform' ]
		else:
			p1 = None
			T0 = meta[ 'transform' ]
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
		if p1 is not None:
			p1 = Proj(p1)
		longs, lats = transform(p1, p1.to_latlong(), eastings, northings)
		return longs, lats
	else:
		BaseException( 'cant reproject to latlong without an input_crs' )

def xyz_to_grid( x, y, z, xi, yi, method='linear', output_dtype=np.float32 ):
	'''
	interpolate points to a grid. simple wrapper around
	matplotlib.mlab.griddata. Points and grid must be
	in the same coordinate system
	x = 1-D np.array of x coordinates / x,y,z must be same length
	y = 1-D np.array of y coordinates / x,y,z must be same length
	z = 1-D np.array of z coordinates / x,y,z must be same length
	grid = tuple of meshgrid as made using `numpy.meshgrid` / `numpy.mgrid`
			order (xi, yi)
	method = 'linear' # hardwired currently due to multiprocessing bug with scipy griddata
	'''
	import numpy as np
	from matplotlib.mlab import griddata
	return griddata( x, y, z, xi, yi, interp=method ).astype( output_dtype )

# def transform_from_latlon( lat, lon ):
# 	''' simple way to make an affine transform from lats and lons coords '''
# 	from affine import Affine
# 	lat = np.asarray( lat )
# 	lon = np.asarray( lon )
# 	trans = Affine.translation(lon[0], lat[0])
# 	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
# 	return trans * scale

# def regrid( x ):
# 	return xyz_to_grid( **x )

def make_df_interpvals( fn ):
	''' this will remove the -9999's from the mix for use in interpolation '''
	from shapely.geometry import Point
	
	# get coords 
	x,y = coordinates( fn )
	# open raster 
	rst = rasterio.open( fn )
	# get data 
	arr = rst.read( 1 )
	
	# make an ndimensional array to leverage numpy indexing speed with x,y,z from the raster 
	new_arr = np.stack([x.ravel(),y.ravel(),arr.ravel()], axis=1)
	new_arr = new_arr[ new_arr[:,2] != -9999,: ] # drop nodata values

	df = pd.DataFrame({'x':new_arr[:,0],'y':new_arr[:,1],'z':new_arr[:,2]})
	df['geometry'] = df.apply( lambda x: Point(x.x,x.y), axis=1 )
	return df


if __name__ == '__main__':
	import os, glob, rasterio
	import numpy as np
	import pandas as pd
	import geopandas as gpd

	# open and set up the alaska data
	ak_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/climatologies/raw/prism/alaska/tmean/ak_tavg_11.txt'
	can_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/climatologies/raw/prism/canada/tmean/caw_tmean_11.asc'
	template_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/templates/akcan_2km/tas_mean_C_ar5_IPSL-CM5A-LR_rcp26_01_2006.tif'
	output_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/climatologies/prism_interpolated'

	# output grid to interpolate to...
	xi, yi = coordinates( template_fn )
	
	# get the AK data prepped to a Pandas DF for interpolation...
	df_ak = make_df_interpvals( ak_fn )
	ak = gpd.GeoDataFrame( df_ak, geometry='geometry', crs={'init':'EPSG:3338'} )
	ak_4326 = ak.to_crs( epsg=4326 )

	# get the AK data prepped to a Pandas DF for interpolation...
	df_can = make_df_interpvals( can_fn )
	can = gpd.GeoDataFrame( df_can, geometry='geometry', crs={'init':'EPSG:4326'} )
	# # # # testing
	# can.to_file( '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/climatologies/test_canada_4326.shp' )
	# # # # 
	can_3338 = can.to_crs( epsg=3338 )
	# scale the canadian data 
	can_3338['z'] = can_3338['z'] / 10

	# concat ak/can_3338 here and prepare for interpolation to the akcan grid
	# df = pd.concat( [ak_4326,can], axis=0 )
	df = pd.concat( [ak,can_3338], axis=0 )
	# update the new values
	df.x = df.geometry.apply( lambda x: x.x )
	df.y = df.geometry.apply( lambda x: x.y )

	# interpolate to the same grid with the np.nan values removed to see if this smoothing will fix it...
	interpolated = xyz_to_grid( np.array(df.x), np.array(df.y), np.array(df.z), xi, yi, method='linear', output_dtype=np.float32 )
	# drop the mask... above function is deprecated in future versions too which stinks...
	interpolated = np.array(interpolated)

	# get some metadata from the template raster for masks and output extent/crs/res/nodata
	with rasterio.open( template_fn ) as tmp:
		meta = tmp.meta.copy()
		meta.update({'compress':'lzw'})
		
		# handle rasterio 1.0 transition from the affine meta keyword.
		if 'affine' in meta.keys():
			_ = meta.pop( 'affine' )

		mask = tmp.read_masks( 1 )
		nodata = tmp.nodata

	# mask it
	interpolated[ mask == 0 ] = nodata

	# write it to disk -- update this to something much better than this shit.
	output_filename = os.path.join(output_dir, os.path.basename(ak_fn).replace('ak', 'akcan')).replace('.txt', '.tif')

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( interpolated , 1 )


