# PRISM CONVERSION FROM ASCII GRIDS -- TASMIN / TASMAX

# header info
# ncols         2015
# nrows         1320
# xllcorner     -2301787.7731349
# yllcorner     108069.7858797
# cellsize      2000
# NODATA_value  -9999

import rasterio, glob, os
from rasterio import Affine
import numpy as np
from pathos import multiprocessing as mp

# input_path = '/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/pr'
# #'/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_2KM_PRISM/Temperature/2km/older'
# output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2'
# groups = ['min_temp', 'max_temp']

# # # STEP 1 -- CONVERT TO GTIFF FROM ASC AND TXT

# list the data we want
variables = [ 'tmin', 'tmax' ]
input_path_ak = '/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_2KM_PRISM/Temperature/2km/older'
input_path_can = '/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_CAN_2km_PRISM/CAN_originals/older'

for variable in variables:
	for ak_test, input_path in zip( [True,False], [input_path_ak,input_path_can] ):
		output_path = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2', variable,'raw_converted' )

		if not os.path.exists( output_path ):
			os.makedirs( output_path )

		if ak_test:
			input_path = input_path_ak
			if variable == 'tmin':
				v = 'min_temp'
			elif variable == 'tmax':
				v = 'max_temp'
			else:
				NotImplemented( 'only tmax / tmin currently supported' )

			files = glob.glob( os.path.join( input_path, v, '*'+variable+'*.txt' ) )
		else:
			input_path = input_path_can
			files = glob.glob( os.path.join( input_path, '*'+variable+'*.asc' ) )

		ext = files[0].split('.')[1]
		output_filenames = [ os.path.join( output_path, os.path.basename( fn ).replace( '.'+ext, '.tif' ) ) for fn in files ]
		crs = {'init':'epsg:4326'}
		args = [ (i,j,crs) for i,j in zip(files, output_filenames) ]

		def bounds_to_extent( bounds ):
			'''
			take input rasterio bounds object and return an extent
			'''
			l,b,r,t = bounds
			return [ (l,b), (r,b), (r,t), (l,t), (l,b) ]
		def convert_to_gtiff( fn, output_filename, crs={'init':'epsg:3338'} ):
			'''
			convert the ascii rasters from PRISM to gtiff
			'''
			print( fn )
			rst = rasterio.open( fn )
			arr = rst.read( 1 ) # get the first and only band
			meta = rst.meta
			meta.update( compress='lzw', driver='GTiff', crs=crs )
			# drop the transform to overcome rasterio warnings
			if 'transform' in meta.keys():
				meta.pop( 'transform' )
			# write them out
			with rasterio.open( output_filename, 'w', **meta ) as out:
				out.write( arr, 1 )
			return output_filename

		if __name__ == '__main__':
			pool = mp.Pool( 32 )
			pool.map( lambda x: convert_to_gtiff( *x ), args )
			pool.close()
			pool.join()

	# # # STEP 2 -- MERGE IT WITH GDAL TOOLS
	# list the data
	caw = sorted( glob.glob( os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2',variable,'raw_converted', 'caw*.tif' ) ) )
	ak = sorted( glob.glob( os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2',variable,'raw_converted', 'ak_*.tif' ) ) )
	grouped = zip( ak, caw )

	# merge these files:
	# log = open( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2/batch_run.bat', 'w' )
	for ak,ca in grouped:
		out = ak.replace( 'ak_', 'akcan_')
		ca_out = ca.replace( '.tif', '_3338.tif' )
		os.system( 'gdalwarp -overwrite -r near -t_srs EPSG:3338 -s_srs EPSG:4326 -ot Float32 ' + ca + ' ' + ca_out )
		ca_scale = ca_out.replace( '.tif', '_scaled.tif' )
		os.system( 'gdal_calc.py --overwrite -A ' + ca_out + ' --outfile=' + ca_scale + ' --calc="A*(0.1)" --NoDataValue=-9999 --type=Float32' )
		os.system( 'gdal_merge.py -init -9999 -n -9999 -a_nodata -9999 -ot Float32 -o ' + out + ' ' + ak + ' ' + ca_scale )
		final = ca.replace( '.tif', '_merged.tif' ).replace( 'raw_converted', 'merged' ).replace( 'caw_', 'akcan_' )
		if not os.path.exists( os.path.dirname(final) ):
			os.makedirs(os.path.dirname(final))
		os.system( 'gdal_translate -co "COMPRESS=LZW" ' + out + ' ' + final )

		# # DUE TO SOME WEIRDNESS WITH VIRTUALENV AND GDAL_MERGE.PY I am writing this out to a text file and running it when not in virtualenv
		# out = ak.replace( 'ak_', 'akcan_')
		# ca_out = ca.replace( '.tif', '_3338.tif' )
		# log.write( 'gdalwarp -overwrite -r near -t_srs EPSG:3338 -s_srs EPSG:4326 -ot Float32 ' + ca + ' ' + ca_out + '\n' )
		# ca_scale = ca_out.replace( '.tif', '_scaled.tif' )
		# log.write( 'gdal_calc.py --overwrite -A ' + ca_out + ' --outfile=' + ca_scale + ' --calc="A*(0.1)" --NoDataValue=-9999 --type=Float32' + '\n' )
		# log.write( 'gdal_merge.py -init -9999 -n -9999 -a_nodata -9999 -ot Float32 -o ' + out + ' ' + ak + ' ' + ca_scale + '\n' )
		# final = ca.replace( '.tif', '_merged.tif' )
		# log.write( 'gdal_translate -co "COMPRESS=LZW" ' + out + ' ' + final + '\n' )


	# # # STEP 3 -- INTERPOLATE / REGRID / MASK to match existing SNAP resources
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
	def xyz_to_grid( x, y, z, grid, method='cubic', output_dtype=np.float32, *args, **kwargs ):
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
		zi = np.flipud( zi ).astype( output_dtype )
		return zi

	import glob, os, rasterio
	import numpy as np
	from rasterio.warp import reproject, RESAMPLING
	import pandas as pd

	# ORIG_AK_RAW = '/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_2KM_PRISM/Temperature/2km/older'

	# TEMPLATE:
	template_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/akcan_template/tas_mean_C_AR5_CCSM4_rcp26_01_2006.tif'
	rst = rasterio.open( template_fn )
	mask = rst.read_masks()

	input_dir = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2', variable, 'merged' )
	files = glob.glob( os.path.join( input_dir, 'akcan_*.tif' ) )

	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	# make an empty raster to hold the new output
	for fn in files:
		template = rasterio.open( template_fn )
		meta = template.meta
		meta.update( compress='lzw' )
		output_filename = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2', variable, os.path.basename( fn ).replace( 'merged', 'final' ) )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( np.empty_like( template.read( 1 ) ), 1 )

		# run it with gdalwarp
		command = 'gdalwarp ' + fn + ' ' + output_filename
		os.system( command )

	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

	# for fn in files:
	# 	cur = rasterio.open( fn )
	# 	# GET THE COORDS AND THE DATA IN A WAY THAT IS INTEPOLATABLE
	# 	lons, lats = coordinates( fn )
	# 	df = pd.DataFrame( {'lons':lons.ravel(), 'lats':lats.ravel(), 'dat':cur.read(1).ravel() } )
	# 	new_df = df[df.dat != -9999]
	# 	new_grid = xyz_to_grid( new_df.lons.tolist(), new_df.lats.tolist(), new_df.dat.tolist(), (lons,lats), method='cubic', output_dtype=np.float32 )
	# 	# new_grid[ np.isnan( new_grid ) ] = -9999
	# 	output_arr = np.empty_like( rst.read( 1 ) )

	# 	# now reproject the new_grid to the extent/resolution/crs of the ALF data -- idiotic crs we use here
	# 	reproject( new_grid, output_arr, src_transform=cur.affine, src_crs={ 'init':'epsg:3338' }, src_nodata=None, \
	# 		dst_transform=rst.affine, dst_crs=rst.crs,\
	# 		dst_nodata=None, resampling=RESAMPLING.cubic_spline, SOURCE_EXTRA=1000 )

	# 	output_arr[ mask == 0 ] = rst.nodata
	# 	new_path = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2', variable, 'prepped' )
	# 	output_filename = os.path.join( new_path, os.path.basename( fn ) )
	# 	meta = rst.meta
	# 	meta.update( compress='lzw' )
	# 	meta.pop( 'transform' )
	# 	with rasterio.open( output_filename, 'w', **meta ) as out:
	# 		out.write( output_arr, 1 )

# # # # END PREP

# THIS IS HOW IT WAS CONVERTED FROM THE TXT/ASC FORMATS
# for group in groups:
# 	# list the data we want to convert to GTiff
# 	# remember that month 14 is the annual average
# 	files = glob.glob( os.path.join( input_path, group, '*.txt' ) )
# 	for fn in files:
# 		print fn
# 		rst = rasterio.open( fn )
# 		arr = rst.read( 1 ) # get the first and only band
# 		output_filename = os.path.join( output_path, os.path.basename(fn).replace( '.txt', '.tif' ) )
# 		meta = rst.meta
# 		meta.update( compress='lzw', driver='GTiff', crs={'init':'epsg:3338'} )
# 		# drop the transform to overcome rasterio warnings
# 		if 'transform' in meta.keys():
# 			meta.pop( 'transform' )
# 		# write them out
# 		with rasterio.open( output_filename, 'w', **meta ) as out:
# 			out.write( arr, 1 )
