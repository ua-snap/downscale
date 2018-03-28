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
	# base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cru_cl20_test_remove'
	# cru_filename = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/grid_10min_pre.dat.gz'
	# variable = 'pre'
	# template_raster_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/templates/akcan_2km/tas_mean_C_ar5_IPSL-CM5A-LR_rcp26_01_2006.tif'
	# # # # # # # # # # # # # #

	# build an output path to store the data generated with this script
	output_path = os.path.join( base_path, 'climatologies','cru_cl20','2km', variable )

	if not os.path.exists( output_path ):
		os.makedirs( output_path )

	months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
	colnames = [ 'lat', 'lon' ] + months
	out_colnames = colnames 

	if variable == 'pre':
		colnames = colnames + [ 'cv{}'.format(m) for m in months ]
	if variable == 'elv':
		colnames = ['lat','lon','01']

	months_lookup = { count+1:month for count, month in enumerate( months ) }
	cru_df = pd.read_csv( cru_filename, delim_whitespace=True, compression='gzip', header=None, names=colnames )
	# slice to the pre values only.  We dont want the cv values for now.
	cru_df = cru_df[ out_colnames ]
	
	# manually flip to PCLL for interpolation
	cru_df['lon'][ cru_df['lon'] < 0 ] = cru_df['lon'][ cru_df['lon'] < 0 ] + 360

	cru_df['geometry'] = cru_df.apply( lambda x: Point( x.lon, x.lat), axis=1 )
	cru_shp = gpd.GeoDataFrame( cru_df, geometry='geometry', crs={'init':'EPSG:4326'} )

	# set bounds to interpolate over
	# xmin, ymin, xmax, ymax = (0,-90, 360, 90)
	xmin, ymin, xmax, ymax = (160, 0, 300, 90) # smaller than global and larger than what we need.
	
	# multiply arcminutes in degree by 360(180) for 10' resolution
	rows = 60 * ( ymax - ymin )
	cols = 60 * ( xmax - xmin )

	# build the output grid
	x = np.linspace( xmin, xmax, cols )
	y = np.linspace( ymax, ymin, rows ) # [ NEW March 2018 ]... flipped ymax/min order in this operation to be north-up...
	xi, yi = np.meshgrid( x, y )

	if variable == 'elv':
		args_list = [{'x':np.array(cru_df['lon']),'y':np.array(cru_df['lat']),'z':np.array(cru_df['01']),'xi':xi,'yi':yi}]
	else:
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

	# TESTING WRITE BELOW
	# with rasterio.open( 'test_pcll_cru_cl20_f.tif','w', **meta ) as out:
	# 	out.write( arr[0,...], 1 )
	# # # # # # # # # # 

	# set up a dir to toss the intermediate files into -- since we are using gdalwarp...
	intermediate_path = os.path.join( output_path, 'intermediates' )
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
	template_meta.update( compress='lzw', nodata=-9999 )
	a,b,c,d = template_raster.bounds
	
	# FLIP IT BACK TO GREENWICH-CENTERED using gdalwarp... then to AKCAN 2km...
	for fn in out_paths:
		os.system( 'gdalwarp -q -co COMPRESS=LZW -overwrite -dstnodata -9999 -multi -t_srs EPSG:4326 -te -180 0 180 90 {} {}'.format( fn, fn.replace( 'PCLL', 'LL' ) ) )
		
		# build an output data set based on the template raster extent and reproject _into_ it
		final_fn = fn.replace( '_PCLL', '' )
		final_fn = os.path.join( output_path, os.path.basename(final_fn) )
		if os.path.exists( final_fn ):
			os.remove( final_fn )

		mask = template_raster.read_masks( 1 ).astype( np.float32 )
		with rasterio.open( final_fn, 'w', **template_meta ) as out:
			out.write( np.empty_like( mask ), 1 )

		os.system( 'gdalwarp -q -wo SOURCE_EXTRA=100 -multi -srcnodata -9999 -dstnodata -9999 {} {}'.format( fn.replace( 'PCLL', 'LL' ), final_fn ) )
		
		# mask newly updated warped output dset
		with rasterio.open( final_fn, 'r+' ) as rst:
			arr = rst.read( 1 )
			arr[ mask == 0 ] = -9999
			
			# round the precip and temperature outputs to the desired precisions
			if variable in ['pre','pr','ppt']:
				arr[ arr != -9999 ] = np.around( arr[ arr != -9999 ], 0 )

			elif variable in ['tmp','tas']:
				arr[ arr != -9999 ] = np.around( arr[ arr != -9999 ], 0 )

			rst.write( arr, 1 )

	print( 'completed run of {}'.format( variable ) )
