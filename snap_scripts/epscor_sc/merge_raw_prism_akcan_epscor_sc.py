# # # # # # # # 
# merge PRISM data from 2016 acquisition from the PRISM Climate Group
# # # # # # # # 

def fill_mask_mismatch( in_arr, value_to_fill=-1000, final_value=-9999, count_missing=(0,1) ):
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
	
	# since we update the array, make sure its not a view...
	arr = np.copy( in_arr ) 
	
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
			vals = np.array([ in_arr[ n ] for n in nlist if in_arr[ n ] != -9999 if in_arr[ n ] != value_to_fill ])
			
			if len( vals ) == 0:
				new_val = final_value
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
	import os, glob, rasterio
	# from rasterio.crs import CRS
	import numpy as np

	base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/raw/prism_raw_2016'

	# get metadata and masks from template raster
	template_rst = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/akcan_template/tas_mean_C_AR5_CCSM4_rcp26_01_2006.tif'
	template = rasterio.open( template_rst )
	xmin, ymin, xmax, ymax = template.bounds
	meta = template.meta
	meta.pop( 'transform' )
	meta.update( compress='lzw', dtype='float32' )
	mask = template.read_masks( 1 ) == 0

	for v in ['tmax', 'tmin']:
		# read files
		ak_files = sorted( glob.glob( os.path.join(base_dir, 'ak', v, '*.asc' ) ) )
		can_files = sorted( glob.glob( os.path.join(base_dir, 'can', v, '*.asc' ) ) )

		akcan_files = zip( ak_files, can_files )

		for ak, can in akcan_files:
			print( 'running: {} - {}'.format( os.path.basename( ak ), os.path.basename( can ) ) )
			output_path = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/merged', v )
			if not os.path.exists( output_path ):
				os.makedirs( output_path )

			# make an extended ak with template extent and fill it with the can data...
			out_fn = os.path.join( output_path, os.path.basename( ak ).replace( '.asc', '.tif' ).replace( 'ak_', 'akcan_') )
			command = 'gdalwarp -multi -q -overwrite -of GTiff -co COMPRESS=LZW -srcnodata -9999 -dstnodata -9999 -tr 2000 2000 -te {} {} {} {} -t_srs EPSG:3338 -s_srs EPSG:3338 {} {}'
			os.system( command.format( xmin, ymin, xmax, ymax, ak, out_fn ) )
			os.system( 'gdalwarp -multi -q -t_srs EPSG:3338 -srcnodata -9999 -s_srs EPSG:4326 {} {}'.format( can, out_fn ) )

			# open the output and read to an array for further filling.
			rst = rasterio.open( out_fn )
			arr = rst.read( 1 )
			rst.close()
			rst = None

			varname_lookup = { 'tmax':'tasmax', 'tmin':'tasmin' }
			out_varname = varname_lookup[ v ]

			final_path = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism', out_varname )
			if not os.path.exists( final_path ):
				os.makedirs( final_path )

			month = os.path.basename( out_fn ).split( '.' )[0].split( '_' )[-1]
			out_name = '_'.join([ out_varname, 'mean', 'C', 'akcan', 'prism', month, '1961', '1990', 'ak83alb' ]) + '.tif'
			output_filename = os.path.join( final_path, out_name )

			# write out filled / masked array
			with rasterio.open( output_filename, 'w', **meta ) as out:
				# print( out.crs )
				# mask the data
				arr[ mask ] = -9999

				# fill missing values in masked extent by giving it `value_to_fill`
				arr[ ( ~mask ) & ( arr == -9999 ) ] = -1000
				filled = fill_mask_mismatch( arr, value_to_fill=-1000, final_value=-9999 )
				filled[ mask ] = -9999
				masked = np.ma.masked_where( filled == -9999, filled, copy=True )
				scaled = np.divide( masked, 10.0 ).astype( np.float32 )
				scaled.fill_value = meta[ 'nodata' ]
				final_arr = np.round( scaled, 1 ).filled()
				out.write( final_arr, 1 )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# TAS / PR FIX Since they have the wrong EXTENT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
import os, glob, rasterio

# get bounds from template
# get metadata and masks from template raster
template_rst = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/akcan_template/tas_mean_C_AR5_CCSM4_rcp26_01_2006.tif'
template = rasterio.open( template_rst )
xmin, ymin, xmax, ymax = template.bounds
meta = template.meta
meta.pop( 'transform' )
meta.update( compress='lzw', dtype='float32' )
mask = template.read_masks( 1 ) == 0

for v in [ 'tas', 'pr' ]:
	data_path = os.path.join( '/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_CAN_2km_PRISM/AK_CAN_geotiffs', v, 'ak83albers' )
	files = sorted( glob.glob( os.path.join( data_path, '*.tif' ) ) )

	output_path = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism', v )
	if not os.path.exists( output_path ):
		os.makedirs( output_path )

	for fn in files:
		print( 'running crop: {}'.format( os.path.basename( fn ) ) )
		out_crop_fn = os.path.join(  '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/cropped', os.path.basename( fn ) )

		if not os.path.exists( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/cropped' ):
			os.makedirs( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/cropped' )

		command = 'gdalwarp -multi -q -overwrite -of GTiff -tr 2000 2000 -dstnodata -3.4e+38 -te {} {} {} {} -co COMPRESS=LZW {} {}'
		os.system( command.format( xmin, ymin, xmax, ymax, fn, out_crop_fn ) )

		# read in the data we just produced
		with rasterio.open( out_crop_fn ) as out:
			arr = out.read( 1 )
		
		# make sure the reference is the same as above.
		out_fn = os.path.join( output_path, os.path.basename( fn ) )
		with rasterio.open( out_fn, 'w', **meta ) as out:
			arr[ mask ] = meta[ 'nodata' ]
			out.write( arr, 1 )

