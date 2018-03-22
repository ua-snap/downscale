# # # PREPROCESS CRU CL20 1961-1990 CLIMATOLOGY DATA (http://www.cru.uea.ac.uk/cru/data/hrg/tmc)
# # author: Michael Lindgren (malindgren@alaska.edu) -- Sept. 2016
# # # #

import numpy as np

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
	# cru_filename = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/grid_10min_sunp.dat.gz'
	# variable = 'sunp'
	# template_raster_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/templates/akcan_10min/akcan_with_nwt_15k_template.tif'

	# # # # # # # # # # # # # #

	# build an output path to store the data generated with this script
	cru_path = os.path.join( base_path, 'climatologies','cru_cl20','10min', variable )

	if not os.path.exists( cru_path ):
		os.makedirs( cru_path )

	months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
	colnames = [ 'lat', 'lon' ] + months
	months_lookup = { count+1:month for count, month in enumerate( months ) }
	cru_df = pd.read_csv( cru_filename, delim_whitespace=True, compression='gzip', header=None, names=colnames )
	
	# manually flip to PCLL for interpolation
	cru_df['lon'][ cru_df['lon'] < 0 ] = cru_df['lon'][ cru_df['lon'] < 0 ] + 360

	cru_df['geometry'] = cru_df.apply( lambda x: Point( x.lon, x.lat), axis=1 )
	cru_shp = gpd.GeoDataFrame( cru_df, geometry='geometry', crs={'init':'EPSG:4326'} )

	# set bounds to interpolate over
	# xmin, ymin, xmax, ymax = (0,-90, 360, 90)
	xmin, ymin, xmax, ymax = (160, 0, 300, 90)

	# multiply arcminutes in degree by 360(180) for 10' resolution
	rows = 60 * ( ymax - ymin )
	cols = 60 * ( xmax - xmin )

	# build the output grid
	x = np.linspace( xmin, xmax, cols )
	y = np.linspace( ymin, ymax, rows )
	xi, yi = np.meshgrid( x, y )
	args_list = [ {'x':np.array(cru_df['lon']),'y':np.array(cru_df['lat']),'z':np.array(cru_df[month]),'xi':xi,'yi':yi} for month in months ]

	# run interpolation in parallel
	interped_grids = mp_map( regrid, args_list, nproc=12 )

	# stack and give a proper nodata value
	arr = np.array([ i.data for i in interped_grids ])
	arr[ np.isnan(arr) ] = -9999
	pcll_affine = transform_from_latlon( y, x )

	meta = {'transform': pcll_affine,
			'count': 1,
			'crs': {'init':'epsg:4326'},
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
		output_filename = os.path.join( intermediate_path, '{}_cru_cl20_akcan_{}_1961-1990_PCLL.tif'.format( variable, months_lookup[ i+1 ] ) )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( arr[ i, ... ], 1 )
		out_paths = out_paths + [ output_filename ]

	# # template dataset
	template_raster = rasterio.open( template_raster_fn )
	resolution = template_raster.res
	template_meta = template_raster.meta
	template_meta.update( compress='lzw', dtype='float32', crs={'init':'epsg:3338'}, nodata=-9999 ) # update some output meta for the new GTiffs
	# if 'transform' in template_meta.keys():
	# 	template_meta.pop( 'transform' )
	a,b,c,d = template_raster.bounds
	
	# FLIP IT BACK TO GREENWICH-CENTERED using gdalwarp... then to AKCAN 10min...
	for fn in out_paths:
		# back to greenwich LL
		os.system( 'gdalwarp -q -overwrite -srcnodata -9999 -dstnodata -9999 -multi -wo SOURCE_EXTRA=100 -t_srs EPSG:4326 -te -180 0 180 90 {} {}'.format( fn, fn.replace( 'PCLL', 'LL' ) ) )
		
		# setup the output filename and pathing
		final_fn = fn.replace( '_PCLL', '' )
		final_fn = os.path.join( cru_path, os.path.basename(final_fn) )
		if os.path.exists( final_fn ):
			os.remove( final_fn )

		# make a new file with the final name that is empty but the same as the template_raster
		mask = template_raster.read_masks( 1 ).astype( np.float32 )
		with rasterio.open( final_fn, 'w', **template_meta ) as out:
			out.write( np.empty_like( mask ), 1 )

		# warp to this new empty dataset
		os.system( 'gdalwarp -q -wo SOURCE_EXTRA=100 -multi -srcnodata -9999 -dstnodata -9999 {} {}'.format( fn.replace( 'PCLL', 'LL' ), final_fn ) )
		with rasterio.open( final_fn, 'r+' ) as rst:
			arr = rst.read( 1 )
			arr[ mask == 0 ] = -9999
			rst.write( arr, 1 )

	print( 'completed run of {}'.format( variable ) )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # EXAMPLE RUN OF THE ABOVE FOR TEM DATA CREATION
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # REGRID CRU-CL20 to SNAP AKCAN 10min
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# import subprocess, os, glob

# script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_10min'
# os.chdir( script_path )
# base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
# cru_filenames = glob.glob('/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/*.dat.gz')
# template_raster_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/templates/akcan_10min/akcan_with_nwt_15k_template.tif'
# for cru_filename in cru_filenames:
# 	print( 'working on: {}'.format( os.path.basename( cru_filename ) ) )
# 	variable = os.path.basename( cru_filename ).split( '_' )[-1].split('.')[0]
# 	done = subprocess.call([ 'ipython', 'cru_cl20_preprocess_10min.py', '--','-p', base_path, '-cru', cru_filename, '-v', variable ,'-tr', template_raster_fn ])