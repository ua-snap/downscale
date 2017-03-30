# # # PREPROCESS CRU CL20 1961-1990 CLIMATOLOGY DATA (http://www.cru.uea.ac.uk/cru/data/hrg/tmc)
# # author: Michael Lindgren (malindgren@alaska.edu) -- March 2017
# # # #

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
			T0 = r.affine  # upper-left pixel corner affine transform
			p1 = Proj( r.crs )
			A = r.read( 1 )  # pixel values

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

def transform_from_latlon( lat, lon ):
	''' simple way to make an affine transform from lats and lons coords '''
	from affine import Affine
	lat = np.asarray( lat )
	lon = np.asarray( lon )
	trans = Affine.translation(lon[0], lat[0])
	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
	return trans * scale

def regrid( x ):
	return xyz_to_grid( **x )

if __name__ == '__main__':
	import os, rasterio
	import numpy as np
	import pandas as pd
	import geopandas as gpd
	from shapely.geometry import Point
	from pathos.mp_map import mp_map
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

	# # # # FOR TESTING # # # #
	# base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
	# cru_filename = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/grid_10min_tmp.dat.gz'
	# variable = 'tmp'
	# template_raster_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/akcan_10min_template/akcan_15k_template.tif'
	# # # # # # # # # # # # #

	# # # THIS IS HARDWIRED CURRENTLY...
	polar_template = rasterio.open( '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/akcan_10min_template/polar_stereo_for_10min_interpolation_extent.tif' )
	#####  ONLY USEFUL FOR 10MIN DATA INCLUDING NWT....

	# build an output path to store the data generated with this script
	cru_path = os.path.join( base_path, 'cru', 'akcan_10min_extent', 'cru_cl20', variable )

	if not os.path.exists( cru_path ):
		os.makedirs( cru_path )

	months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
	colnames = [ 'lat', 'lon' ] + months
	months_lookup = { count+1:month for count, month in enumerate( months ) }
	cru_df = pd.read_csv( cru_filename, delim_whitespace=True, compression='gzip', header=None, names=colnames )
	
	# set bounds to interpolate over
	# xmin, ymin, xmax, ymax = (0,-90, 360, 90)
	xmin, ymin, xmax, ymax = (0, 20, 360, 90)

	# drop everything below the ymin...
	# manually flip to PCLL for interpolation
	cru_df = cru_df[ cru_df.lat >= ymin ]

	# convert to a spatial data frame
	xy = zip(cru_df['lon'].tolist(), cru_df['lat'].tolist())
	points = mp_map( lambda i: Point(i), xy, nproc=50 )
	cru_df[ 'geometry' ] = points
	cru_shp = gpd.GeoDataFrame( cru_df, geometry='geometry', crs={'init':'EPSG:4326'} )

	# polar coords for interpolation...
	p_xmin, p_ymin, p_xmax, p_ymax = list( polar_template.bounds )

	# flip it to polar stereo for bettter (not perfect) interpolation
	cru_shp_polar = cru_shp.to_crs( epsg=3413 )

	# HARDWIRED FROM TEMPLATE EXTENT...
	rows, cols = polar_template.shape
	
	# build the output grid
	xi, yi = coordinates( polar_template.name )

	# build args to run in parallel
	lons = [ geom.x for geom in cru_shp_polar.geometry ]
	lats = [ geom.y for geom in cru_shp_polar.geometry ]
	args_list = [ {'x':np.array(lons),'y':np.array(lats),'z':np.array(cru_shp_polar[month]),'xi':xi,'yi':yi} for month in months ]

	# run interpolation in parallel
	interped_grids = mp_map( regrid, args_list, nproc=12 )

	# stack and give a proper nodata value
	arr = np.array([ i.data for i in interped_grids ])
	arr[ np.isnan(arr) ] = -9999

	meta = {'affine': polar_template.affine,
			'count': 1,
			'crs': {'init':'epsg:3413'},
			'driver': u'GTiff',
			'dtype': 'float32',
			'height': rows,
			'nodata': -9999,
			'width': cols,
			'compress':'lzw'}

	# set up a dir to toss the intermediate files into -- since we are using gdalwarp...
	intermediate_path = os.path.join( cru_path, 'intermediates' )
	if not os.path.exists( intermediate_path ):
		os.makedirs( intermediate_path )

	out_paths = []
	for i in range( arr.shape[0] ):
		output_filename = os.path.join( intermediate_path, '{}_cru_cl20_akcan_{}_1961-1990_POLARSTEREO.tif'.format( variable, months_lookup[ i+1 ] ) )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( arr[ i, ... ], 1 )
		out_paths = out_paths + [ output_filename ]

	# template dataset for warping into 
	template_raster = rasterio.open( template_raster_fn )
	resolution = template_raster.res
	template_meta = template_raster.meta
	template_meta.update( compress='lzw', dtype='float32', crs={'init':'epsg:3338'}, nodata=-9999 ) # update some output meta for the new GTiffs
	if 'transform' in template_meta.keys():
		template_meta.pop( 'transform' )
	a,b,c,d = template_raster.bounds
	template_arr = template_raster.read( 1 ).astype( np.float32 )

	# now lets try to reproject it to 3338 using the template raster...
	for fn in out_paths:
		# FIRST REPROJECT THE 3413 DATA TO 4326 @~10'
		print( fn )
		out_fn = fn.replace( '_POLARSTEREO.tif', '_GCLL.tif' )
		command = 'gdalwarp -srcnodata -9999 -dstnodata -9999 -overwrite -wo SOURCE_EXTRA=10000 -wo NUM_THREADS=5 -multi -r cubic -s_srs EPSG:3413 -t_srs EPSG:4326 -tr 0.16667 0.16667 -te -180 20 180 90 {} {}'.format( fn, out_fn )
		os.system( command )
		
		# # flip to PCLL? 
		# out_fn_pcll = out_fn.replace('GCLL', 'PCLL')
		# command = 'gdalwarp -multi -overwrite -wo SOURCE_EXTRA=1000 -wo NUM_THREADS=5 -t_srs WGS84 -te 0 20 360 90 {} {}'.format( out_fn, out_fn_pcll )
		# os.system( command )

		# make an empty layer with the same meta as template raster for warping into.
		output_filename = os.path.join( cru_path, os.path.basename( out_fn.replace( '_GCLL.tif', '.tif' ) ) )
		with rasterio.open( output_filename, 'w', **template_meta ) as rst:
			rst.write( np.empty_like( template_arr ), 1 )

		# # reproject TO this new empty layer... EPSG:3338
		command = 'gdalwarp -r cubic -wo SOURCE_EXTRA=1000 -wo NUM_THREADS=5 -multi -srcnodata -9999 {} {}'.format( out_fn, output_filename )
		os.system( command )

		# mask it
		with rasterio.open( output_filename ) as rst:
			out_arr = rst.read( 1 )
			out_arr[ template_arr == 0 ] = -9999

		with rasterio.open( output_filename, 'w', **template_meta ) as out:
			out.write( out_arr, 1 )
		break

	print( 'completed run of {}'.format( variable ) )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # EXAMPLE RUN OF THE ABOVE FOR TEM DATA CREATION
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # REGRID CRU-CL20 to SNAP AKCAN 2KM
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# import subprocess, os, glob

# script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_10min'
# os.chdir( script_path )
# base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
# cru_filenames = [ fn for fn in glob.glob('/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/*.dat.gz') if 'tmp' in fn or 'pre' in fn ]
# template_raster_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/akcan_10min_template/akcan_15k_template.tif'
# for cru_filename in cru_filenames:
# 	print( 'working on: {}'.format( os.path.basename( cru_filename ) ) )
# 	variable = os.path.basename( cru_filename ).split( '_' )[-1].split('.')[0]
# 	done = subprocess.call([ 'ipython', 'cru_cl20_preprocess.py', '--','-p', base_path, '-cru', cru_filename, '-v', variable ,'-tr', template_raster_fn ])
