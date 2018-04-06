# # # PREPROCESS CRU CL20 1961-1990 CLIMATOLOGY DATA (http://www.cru.uea.ac.uk/cru/data/hrg/tmc)
# # author: Michael Lindgren (malindgren@alaska.edu) -- Sept. 2016
# # # #

import numpy as np
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
	import os, rasterio, glob
	import numpy as np
	import pandas as pd
	import geopandas as gpd
	from shapely.geometry import Point
	from pathos.mp_map import mp_map
	from downscale import utils
	import argparse

	# # parse the commandline arguments
	# parser = argparse.ArgumentParser( description='preprocess WorldClim 1.4 data to the AKCAN extent required by SNAP' )
	# parser.add_argument( "-p", "--base_path", action='store', dest='base_path', type=str, help="path to parent directory with a subdirector(ies)y storing the data" )
	# parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="string abbreviated name of the variable being processed." )
	# parser.add_argument( "-tr", "--template_raster_fn", action='store', dest='template_raster_fn', type=str, help="string path to a template raster dataset to match the WorldClim 1.4 to..." )

	# # parse and unpack the args
	# args = parser.parse_args()
	# base_path = args.base_path
	# variable = args.variable
	# template_raster_fn = args.template_raster_fn

	# # # FOR TESTING # # # #
	base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/worldclim'
	variable = 'tmean'
	template_raster_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/akcan_template/tas_mean_C_AR5_CCSM4_rcp26_01_2006.tif'
	# # # # # # # # # # # # #
	
	# get the files and sort them chronologically
	files = glob.glob( os.path.join( base_path, 'raw', variable, '*.bil' ) )
	files = pd.Series( { int(fn.split('.')[0].split('_')[-1]):fn for fn in files } )
	files = files.sort_index().tolist()

	# build an output path to store the data generated with this script
	wc_path = os.path.join( base_path, 'worldclim_1_4', variable )

	if not os.path.exists( wc_path ):
		os.makedirs( wc_path )

	# DO THE RIDICULOUSLY SIMPLY THING OF REGRIDDING THE DATA FROM GLOBAL LL GREENWICH TO 3338 2KM MASKED
	template_raster = rasterio.open( template_raster_fn )
	meta = template_raster.meta
	meta.update( compress='lzw' )
	if 'transform' in meta.keys():
		_ = meta.pop( 'transform' )

	for fn in files:
		print( fn )
		with rasterio.open( fn ) as rst:
			rst_meta = rst.meta
			rst_meta.update( compress='lzw', dtype='float32' )

			arr = rst.read( 1 ).astype( np.float32 )
			arr[ arr != -9999 ] = arr[ arr != -9999 ] * .10
			# new name for scaled values raster
			variable, month = os.path.basename( fn ).split( '.' )[0].split( '_' )

			intermediate_path = os.path.join( wc_path, 'intermediates' )
			if not os.path.exists( intermediate_path ):
				os.makedirs( intermediate_path )

			new_fn = os.path.join( intermediate_path, '{}_{}_worldclim_1_4_{}_1961-1990.tif'.format( variable, 'mean', month ) )
			with rasterio.open( new_fn, 'w', **rst_meta ) as out:
				out.write( arr, 1 )

		output_filename = os.path.join( wc_path, os.path.basename( new_fn ).replace( '.tif', '_akcan.tif' ) )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( np.empty_like( template_raster.read( 1 ) ), 1 )

		command = 'gdalwarp -multi -r bilinear -s_srs EPSG:4326 -wo SOURCE_EXTRA=100 {} {}'.format( new_fn, output_filename )
		os.system( command )

		with rasterio.open( output_filename, mode='r+') as out:
			arr = out.read( 1 )
			mask = template_raster.read_masks(1)
			arr[ mask == 0 ] = template_raster.nodata
			out.write( arr, 1 )

# # # # # # END FOR THE EASY SOLUTION...
# OLDER BELOW, BUT GOOD STUFF

# 	# get centroid coordinates as a meshgrid
# 	# bounds = (160, 0, 300, 90)
# 	bounds = (-180, 0, 180, 90)
# 	rst = rasterio.open( files[0] )
# 	window = rst.window( *bounds )
# 	window_affine = rst.window_transform( window )

	
# 	# 2. Rotate to PCLL and clip to a decent extent for processing.
# 	# # # NOTE ROTATE TO PCLL USING GDALWARP FIRST, THEN DO THE OTHER STUFF.
# 	arr = np.array( [ rasterio.open( fn ).read( 1, window=window ) for fn in files ] )
# 	# # # # WHAT IF INSTEAD OF DOING THIS PROPER COORDINATE CREATION FOR THE INTERPOLATION WHICH IS ESSENTIALLY MEANINGLESS
# 	# WHY NOT CREATE A NEW SIMPLE COORDINATES WITH np.linspace and the nrows and ncols.  Then we just put it back in the right way. 
# 	# THE CRS IS NOT REALLY BEING TAKEN INTO ACCOUNT IN THE INTERPOLATION ANYWAY.
# 	lons, lats = coordinates( meta={'affine':window_affine}, numpy_array=rst.read(1, window=window), input_crs={'init':'epsg:4326'} )
# 	arr_pcll, lons = utils.shiftgrid( 0., arr, lons )

# 	# 	- maybe convert to points here and flip the coords manually in a DataFrame?
# 	months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
# 	# d.update( lat=lat, lon=lon )
# 	# wc_df = pd.DataFrame( d )
# 	# # 	- drop the rows with nodata values:
# 	# wc_df. # << - this is the most important peice to make this run.

# 	# # 	- rotate that grid to pcll
# 	# wc_df['lon'][ wc_df['lon'] < 0 ] = wc_df['lon'][ wc_df['lon'] < 0 ] + 360

# 	# THESE ARE POTENTIALLY DONE BELOW
# 	# 3. interpolate across space (xyz_to_grid)
# 	# 4. Rotate back to Greenwich
# 	# 5. regrid to AKCAN 2km
# 	# 6. mask and write out to GTiff.
	
# 	months_lookup = { count+1:month for count, month in enumerate( months ) }

# 	# cru_df['geometry'] = cru_df.apply( lambda x: Point( x.lon, x.lat), axis=1 )
# 	# cru_shp = gpd.GeoDataFrame( cru_df, geometry='geometry', crs={'init':'EPSG:4326'} )

# 	# set bounds to interpolate over
# 	# xmin, ymin, xmax, ymax = (0,-90, 360, 90)
# 	xmin, ymin, xmax, ymax = bounds

# 	# multiply arcminutes in degree by 360(180) for 10' resolution
# 	rows = 60 * ( ymax - ymin )
# 	cols = 60 * ( xmax - xmin )

# 	# build the output grid
# 	x = np.linspace( xmin, xmax, cols )
# 	y = np.linspace( ymin, ymax, rows )
# 	xi, yi = np.meshgrid( x, y )
# 	args_list = [ {'x':np.array(wc_df['lon']),'y':np.array(wc_df['lat']),'z':np.array(wc_df[month]),'xi':xi,'yi':yi} for month in months ]

# 	# run interpolation in parallel
# 	interped_grids = mp_map( regrid, args_list, nproc=12 )

# 	# stack and give a proper nodata value
# 	arr = np.array([ i.data for i in interped_grids ])
# 	arr[ np.isnan(arr) ] = -3.4e+38
# 	pcll_affine = transform_from_latlon( y, x )

# 	meta = {'affine': pcll_affine,
# 			'count': 1,
# 			'crs': {'init':'epsg:4326'},
# 			'driver': u'GTiff',
# 			'dtype': 'float32',
# 			'height': rows,
# 			'nodata': -3.4e+38,
# 			'width': cols,
# 			'compress':'lzw'}

# 	# set up a dir to toss the intermediate files into -- since we are using gdalwarp...
# 	intermediate_path = os.path.join( cru_path, 'intermediates' )
# 	if not os.path.exists( intermediate_path ):
# 		os.makedirs( intermediate_path )

# 	out_paths = []
# 	for i in range( arr.shape[0] ):
# 		output_filename = os.path.join( intermediate_path, '{}_cru_cl20_akcan_{}_1961-1990_PCLL.tif'.format( variable, months_lookup[ i+1 ] ) )
# 		with rasterio.open( output_filename, 'w', **meta ) as out:
# 			out.write( arr[ i, ... ], 1 )
# 		out_paths = out_paths + [ output_filename ]

# 	# # template dataset
# 	template_raster = rasterio.open( template_raster_fn )
# 	resolution = template_raster.res
# 	template_meta = template_raster.meta
# 	template_meta.update( compress='lzw' )
# 	a,b,c,d = template_raster.bounds
	
# 	# FLIP IT BACK TO GREENWICH-CENTERED using gdalwarp... then to AKCAN 2km...
# 	for fn in out_paths:
# 		os.system( 'gdalwarp -overwrite -dstnodata -3.4e+38 -multi -t_srs EPSG:4326 -te -180 0 180 90 {} {}'.format( fn, fn.replace( 'PCLL', 'LL' ) ) )
		
# 		final_fn = fn.replace( '_PCLL', '' )
# 		final_fn = os.path.join( cru_path, os.path.basename(final_fn) )
# 		if os.path.exists( final_fn ):
# 			os.remove( final_fn )

# 		mask = template_raster.read_masks( 1 ).astype( np.float32 )
# 		with rasterio.open( final_fn, 'w', **template_meta ) as out:
# 			out.write( np.empty_like( mask ), 1 )

# 		os.system( 'gdalwarp -wo SOURCE_EXTRA=100 -multi -srcnodata -3.4e+38 -dstnodata -3.4e+38 {} {}'.format( fn.replace( 'PCLL', 'LL' ), final_fn ) )
# 		# os.system( 'gdalwarp -overwrite -t_srs EPSG:3338 -co COMPRESS=LZW -wo SOURCE_EXTRA=100 -multi -srcnodata {} -dstnodata {} {} {}'.format( -3.4e+38, -3.4e+38, fn.replace( 'PCLL', 'LL' ), final_fn ) )
# 		with rasterio.open( final_fn, 'r+' ) as rst:
# 			arr = rst.read( 1 )
# 			arr[ mask == 0 ] = -3.4e+38
# 			rst.write( arr, 1 )

# 	print( 'completed run of {}'.format( variable ) )



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # EXAMPLE RUN OF THE ABOVE FOR TEM DATA CREATION
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # REGRID CRU-CL20 to SNAP AKCAN 2KM
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# import subprocess, os, glob

# script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/tem_inputs_iem'
# os.chdir( script_path )
# base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/cru'
# cru_filenames = glob.glob('/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/*.dat.gz')
# template_raster_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/akcan_template/tas_mean_C_AR5_CCSM4_rcp26_01_2006.tif'
# for cru_filename in cru_filenames:
# 	print( 'working on: {}'.format( os.path.basename( cru_filename ) ) )
# 	variable = os.path.basename( cru_filename ).split( '_' )[-1].split('.')[0]
# 	done = subprocess.call([ 'ipython', 'cru_cl20_preprocess.py', '--','-p', base_path, '-cru', cru_filename, '-v', variable ,'-tr', template_raster_fn ])



