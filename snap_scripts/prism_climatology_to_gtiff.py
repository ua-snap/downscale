# # PRISM CONVERSION TO SNAP -- TMIN / TMAX 

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# # NOTES:
# --------
# 	header info
# 	ncols         2015
# 	nrows         1320
# 	xllcorner     -2301787.7731349
# 	yllcorner     108069.7858797
# 	cellsize      2000
# 	NODATA_value  -9999
#
# # DATA GOTCHAS:
# -------------
# 	The SNAP raw data holdings for prism climatologies are in 2 groups W.Canada & AK
#   # ORIGINALLY WE WORKED WITH THESE: but we have since switched to the 2016 version
# 	- [AK] '/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_2KM_PRISM/Temperature/2km/older'
#	- [CAW] '/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_CAN_2km_PRISM/CAN_originals/older'
#
#	To further complicate things, the AK raw data is in 3338 and the CAW data is in 4326. and the months are identified
#	differently with either 1 digit months or 2 digit months
#	
# 	There is also a strange data point in the 2016 version of the data that has values of -100 on the map in very odd locations
#	like SEAK in January, which doesnt make any sense.  I am going to have to fill these in with the average of the neighbors 
# 	since they are directly on the land-sea divide.  This is a very odd thing and I will inform TK about it.
#
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

import numpy as np

def bounds_to_extent( bounds ):
	'''
	take input rasterio bounds object and return an extent
	'''
	l,b,r,t = bounds
	return [ (l,b), (r,b), (r,t), (l,t), (l,b) ]
def convert_to_gtiff( fn, output_filename, input_crs={'init':'epsg:3338'} ):
	'''
	convert the ascii rasters from PRISM to gtiff
	'''
	import subprocess
	output_crs = '+datum=NAD83 +lat_0=50 +lat_1=55 +lat_2=65 +lon_0=-154 +no_defs +proj=aea +units=m +x_0=0 +y_0=0' # 'EPSG:3338'
	subprocess.call([ 'gdalwarp', '-q', '-overwrite' ,'-of', 'GTiff', '-tap', '-tr', '2000', '2000', '-s_srs', \
						input_crs.values()[0].upper(), '-srcnodata', '-9999', '-dstnodata', '-9999', \
						'-t_srs', output_crs, fn, output_filename ])
	return output_filename
def convert_fill( fn, output_filename, input_crs={'init':'epsg:3338'}, value_to_fill=-100, final_value=-9999, count_missing=(0,1)):
	fn = convert_to_gtiff( fn, output_filename, input_crs )
	with rasterio.open( fn, 'r+' ) as out:
		in_arr = out.read( 1 )
		out.write( fill_mask_mismatch( in_arr, value_to_fill, final_value ), 1 )
	return output_filename
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
	# zi = np.flipud( zi ).astype( output_dtype )
	return zi.astype( output_dtype )
def fill_mask_mismatch( in_arr, value_to_fill=-100, final_value=-9999, count_missing=(0,1) ):
	'''
	recurse through the sea ice array and the error locations array
	and fill in the missing data that is not inland as the mean of 
	the surrounding (queens case) pixels truncated to an integer.

	dont mind the count_missing variable, it is simply a placeholder 
	that is used in storing some information during the recursion. 
	I'm sure there is better practice here but this is elegant enough
	for this application (small raster sizes).

	all other values inland are returned as value 255

	returns a filled sea ice array that can be easily written to disk
	using rasterio.

	'''
	height, width = in_arr.shape
	arr = np.copy( in_arr ) # since we update the array, make sure its not a view...
	# grab the index of the error locations
	ind = np.where( arr == value_to_fill )
	ind = zip( *ind )
	count1, count2 = count_missing
	# print count_missing
	if count1 == count2:
		# fill in the remainders with 255 and return
		for ii in ind:
			in_arr[ ii ] = final_value
		return in_arr
	else:
		# fill em if we can, if not, make em 255 -- we could use more neighbors but that sounds hairy
		# all_missing_neighbors = {(i,j):[ (i-1,j+0), (i+0,j-1), (i+0,j+1), (i+1,j+0) ] for i,j in ind } # rooks
		all_missing_neighbors = { (i,j):[ (i-1,j+0), (i+0,j-1), (i+0,j+1), (i+1,j+0), (i+1,j+1), \
									(i-1,j+1), (i-1,j-1), (i+1,j-1) ] for i,j in ind } # queens

		for missing, neighbors_list in all_missing_neighbors.iteritems():
			nlist = neighbors_list # just hold the name for testing
			nlist = [ (i,j) for i,j in nlist if i >= 0 and i < height if j >= 0 and j < width ]
			vals = np.array([ in_arr[ n ].tolist() for n in nlist if n != -9999 if in_arr[ n ] <= -1000 ]) #if sic_arr[ n ] != 0 
				
			vals = vals[ vals >= 0 ].tolist()
			if len( vals ) == 0:
				new_val = -9999
				arr[ missing ] = 1 # keep flag for missing
			elif len( vals ) > 0 :
				new_val = int( np.mean( vals ) )
				arr[ missing ] = 0 # remove the flag for missing
				in_arr[ missing ] = new_val
			else:
				print 'error'

		count_missing = ( count_missing[1], len( arr[ arr == 1 ] ) )
		return fill_mask_mismatch( in_arr, value_to_fill, final_value, count_missing )

if __name__ == '__main__':
	import rasterio, glob, os
	from rasterio import Affine
	from rasterio.warp import reproject, RESAMPLING
	import numpy as np
	from pathos import multiprocessing as mp
	import pandas as pd

	# list the data we want
	variables = [ 'tmin', 'tmax' ]
	regions = ['ak', 'can']
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data'
	template_fn = os.path.join( base_path, 'akcan_template', 'tas_mean_C_AR5_CCSM4_rcp26_01_2006.tif' )
	template = rasterio.open( template_fn )

	# # # # # # # # # # # # # # # # # # # # # # # # # #
	# # # STEP 1 -- CONVERT TO GTIFF FROM ASC AND TXT #
	# # # # # # # # # # # # # # # # # # # # # # # # # #
	print( 'STEP 1' )
	for variable in variables:
		for region in regions:
			raw_path = os.path.join( base_path, 'prism_raw_2016', region, variable )
			output_path = os.path.join( base_path, 'prism_v2', variable, 'raw_converted' )

			if not os.path.exists( output_path ):
				os.makedirs( output_path )

			# different raw crs's
			if region == 'ak':
				crs = { 'init':'epsg:3338' }
			elif region == 'can':
				crs = { 'init':'epsg:4326' }
			else:
				NotImplemented( 'only supported regions are `ak` or `can`' )

			files = glob.glob( os.path.join( raw_path, '*.asc' ) )
			ext = files[0].split('.')[1]
			output_filenames = [ os.path.join( output_path, os.path.basename( fn ).replace( '.'+ext, '.tif' ) ) for fn in files ]
			value_to_fill = -1000
			final_value = -9999
			args = [ (i,j,crs,value_to_fill,final_value) for i,j in zip(files, output_filenames) ]

			# run in parallel	
			pool = mp.Pool( 32 )
			# pool.map( lambda x: convert_to_gtiff( *x ), args )
			pool.map( lambda x: convert_fill( *x ), args )
			pool.close()
			pool.join()


		# # # # # # # # # # # # # # # # # # # # # # # 
		# # # STEP 2 -- MERGE IT WITH GDAL TOOLS
		# # # # # # # # # # # # # # # # # # # # # # # 
		print( 'STEP 2' )
		# list the converted files
		can = sorted( glob.glob( os.path.join( output_path, 'caw*.tif' ) ) )
		ak = sorted( glob.glob( os.path.join( output_path, 'ak_*.tif' ) ) )
		grouped = zip( ak, can )

		# merge these files:
		for ak,can in grouped:
			out = ak.replace( 'ak_', 'akcan_')
			ak_scale = ak.replace( '.tif', '_scaled.tif' )
			can_scale = can.replace( '.tif', '_scaled.tif' )
			os.system( 'gdal_calc.py --overwrite -A ' + ak + ' --outfile=' + ak_scale + ' --calc="A/10" --NoDataValue=-9999 --type=Float32' ) # --creation-option="COMPRESS=LZW"
			os.system( 'gdal_calc.py --overwrite -A ' + can + ' --outfile=' + can_scale + ' --calc="A/10" --NoDataValue=-9999 --type=Float32' ) # --creation-option="COMPRESS=LZW"
			
			# get template raster bounds to pass to gdal_merge.py
			l, b, r, t = template.bounds # ulx uly lrx lry

			# merge it
			os.system( 'gdal_merge.py -q -ps 2000 2000 -tap -init -9999 -n -9999 -a_nodata -9999 -ot Float32 -ul_lr %s %s %s %s -o %s %s %s' % (l,t,r,b,out,ak_scale,can_scale) )

			# # fill in the -100 data that is in error along the coastlines -- 2016 version
			# with rasterio.open( out, 'r+' ) as rst:
			# 	rst.write( fill_mask_mismatch( rst.read( 1 ), -100, -9999 ), 1 )

			final = can.replace( '.tif', '_merged.tif' ).replace( 'raw_converted', 'merged' ).replace( 'caw_', 'akcan_' )
			if not os.path.exists( os.path.dirname( final ) ):
				os.makedirs( os.path.dirname( final ) )
			
			os.system( 'gdal_translate -q -co "COMPRESS=LZW" ' + out + ' ' + final )
		
		# # # # # # # # # # # # # # # # # # # # # # # 
		# # # STEP 3 -- INTERPOLATE / REGRID / MASK to 
		#				match existing SNAP resources
		# # # # # # # # # # # # # # # # # # # # # # # 
		print( 'STEP 3' )
		mask = template.read_masks( 1 )
		input_dir = os.path.join( base_path, 'prism_v2', variable, 'merged' )
		files = glob.glob( os.path.join( input_dir, 'akcan_*.tif' ) )
		lons, lats = coordinates( files[0] )
		t_lons, t_lats = coordinates( template_fn )
		
		for fn in files:
			print( fn )
			cur = rasterio.open( fn )
			df = pd.DataFrame( {'lons':lons.ravel(), 'lats':lats.ravel(), 'dat':cur.read(1).ravel() } )
			new_df = df[ df.dat != -9999 ]
			new_grid = xyz_to_grid( new_df.lons.tolist(), new_df.lats.tolist(), new_df.dat.tolist(), (t_lons,t_lats), method='cubic', output_dtype=np.float32 )
			new_grid[ (mask == 0) ] = template.nodata
			new_grid[ np.isnan( new_grid ) ] = -9999
			
			# output_arr = np.empty_like( template.read( 1 ) )

			# # # now reproject the new_grid to the extent/resolution/crs of the ALF data -- idiotic crs we use here
			# reproject( new_grid, output_arr, src_transform=cur.affine, src_crs={ 'init':'epsg:3338' }, src_nodata=-9999, \
			# 			dst_transform=template.affine, dst_crs=template.crs,\
			# 			dst_nodata=template.nodata, resampling=RESAMPLING.cubic_spline, SOURCE_EXTRA=1000 )

			# output_arr[ (mask == 0) ] = template.nodata
			# output_arr[ (mask == 0) | (output_arr == -9999) ] = template.nodata

			# # fill in the -9999 data that is in error along the coastlines -- 2016 version
			# new_grid = fill_mask_mismatch( new_grid, -9999, template.nodata )
			output_arr = new_grid
			
			new_path = os.path.join( base_path, 'prism_v2', variable )
			if not os.path.exists( new_path ):
				os.makedirs( new_path )

			output_filename = os.path.join( new_path, os.path.basename( fn ) )
			meta = template.meta
			meta.update( compress='lzw' )
			meta.pop( 'transform' )
			with rasterio.open( output_filename, 'w', **meta ) as out:
				out.write( output_arr, 1 )

			with rasterio.open( output_filename, 'r+' ) as out:
				mask2 = out.read_masks( 1 )
				arr = out.read( 1 )
				arr[ (mask != 0) & (arr == -3.39999995e+38) ] = -9999
				arr = fill_mask_mismatch( arr, -9999, template.nodata )
				out.write( arr, 1 )
