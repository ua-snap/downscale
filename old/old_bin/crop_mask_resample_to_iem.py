def resample_to_1km( x, template_raster_mask ):
	'''
	template_raster_mask should be a mask in in the res/extent/origin/crs of the 
	existing TEM IEM products.
	'''
	import rasterio, os
	from rasterio.warp import RESAMPLING, reproject
	import numpy as np

	fn = os.path.basename( x )
	fn_split = fn.split( '.' )[0].split( '_' ) 
	if '_cru_' in fn:
		output_path = os.path.dirname( x ).replace( '/cru_ts31/', '/IEM/cru_ts31/' ) # hardwired!
		fn_parts = ['variable', 'metric', 'model_1', 'model_2', 'kind', 'month', 'year']
		fn_dict = dict( zip( fn_parts, fn_split ) )
		fn_dict.update( scenario='historical', model='cru_ts31' )
	else:
		output_path = os.path.dirname( x ).replace( '/ar5/', '/IEM/ar5/' ) # hardwired!
		fn_parts = ['variable', 'metric', 'model', 'scenario', 'ensemble', 'month', 'year']
		fn_dict = dict( zip( fn_parts, fn_split ) )

	try:
		if not os.path.exists( output_path ):
			os.makedirs( output_path )
	except:
		pass

	fn_switch = { 'cld':'_'.join([ 'cld','mean','pct','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif',
		'vap':'_'.join(['vap','mean','hPa','iem', fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif', 
		'tas':'_'.join(['tas','mean','C','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif',
		'hur':'_'.join(['hur','mean','pct','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif' }

	output_filename = os.path.join( output_path, fn_switch[ fn_dict[ 'variable' ] ] )

	rst = rasterio.open( x )
	rst_arr = rst.read( 1 )

	template_arr = template_raster_mask.read( 1 )
	template_meta = template_raster_mask.meta
	template_meta.update( compress='lzw', nodata=rst.nodata )
	
	if 'transform' in template_meta.keys():
		template_meta.pop( 'transform' )

	output_arr = np.empty_like( template_arr.astype( np.float32 ) )
	output_arr[ template_arr == 0 ] = rst.nodata

	src_crs = {'init':'epsg:3338'}
	dst_crs = {'init':'epsg:3338'}

	reproject( rst_arr, output_arr, src_transform=rst.affine, src_crs=src_crs, src_nodata=rst.nodata, \
			dst_transform=template_raster_mask.affine, dst_crs=dst_crs,\
			dst_nodata=rst.nodata, resampling=RESAMPLING.cubic_spline, num_threads=2 )

	with rasterio.open( output_filename, 'w', **template_meta ) as out:
		output_arr[ template_arr == 0 ] = rst.nodata
		out.write( output_arr, 1 )
	return output_filename


if __name__ == '__main__':
	import os, glob, rasterio
	import numpy as np
	import pandas as pd
	from functools import partial
	from pathos import multiprocessing as mp

	# some setup:
	input_path = '/Data/malindgren/cru_november_final/ar5' # '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_november_final/ar5'
	template_raster_mask_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/extents/IEM_Mask_1km.tif'
	ncores = 20

	# read in the template raster mask must be 0-nodata, 1-data
	template_raster_mask = rasterio.open( template_raster_mask_fn )
	resample_to_1km_partial = partial( resample_to_1km, template_raster_mask=template_raster_mask )
	
	models = [ 'IPSL-CM5A-LR', 'GISS-E2-R', 'MRI-CGCM3', 'CCSM4', 'GFDL-CM3' ]
	variables = ['cld', 'vap' ] # ['tas'] # ['hur'] # run the HUR after the cld/vap which is required for TEM
	# tas_files = sorted( glob.glob( os.path.join( tas_input_path, model, 'tas', 'downscaled', '*.tif' ) ) )

	path_list = [ os.path.join( input_path, model, variable, 'downscaled', '*.tif' ) for model in models for variable in variables ]

	# temporary for failes EOS run:
	complete = '/Data/malindgren/cru_november_final/ar5/IPSL-CM5A-LR/cld/downscaled/*.tif'
	path_list = [ path for path in path_list if path != complete ]
	# end temporary

	for path in path_list:
		print path 
		files = glob.glob( path )

		# run it in parallel
		pool = mp.Pool( processes=ncores )
		pool.map( lambda x: resample_to_1km_partial( x=x ), files )
		pool.close()



	# for root, subs, files in os.walk( input_path ):
	# 	if root.endswith( 'downscaled' ):
	# 		print 'running: %s' % root
	# 		# add back in the file paths from the root
	# 		files = [ os.path.join( root, i ) for i in files ]

	# 		# run it in parallel
	# 		pool = mp.Pool( processes=ncores )
	# 		pool.map( lambda x: resample_to_1km_partial( x=x ), files )
	# 		pool.close()



# # # # MAKE A MASK TO USE AS THE TEMPLATE RASTER
# template_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/extents/tas_mean_C_iem_cru_TS31_01_1901.tif'
# template = rasterio.open( template_fn )
# template_meta = template.meta
# template_mask = template.read_masks( 1 )
# template_arr = template.read( 1 )
# template_arr[ template_mask != 0 ] = 1
# template_arr[ template_mask == 0 ] = 0

# template_meta.update( compress='lzw', crs={'init':'epsg:3338'} )
# output_filename = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/extents/IEM_Mask_1km.tif'
# with rasterio.open( output_filename, 'w', **template_meta ) as out:
# 	out.write( template_arr, 1 )

# # # # # 
# def standardized_fn_to_vars( fn ):
# 	''' take a filename string following the convention for this downscaling and break into parts and return a dict'''
# 	fn = os.path.basename( fn )
# 	fn_split = fn.split( '.' )[0].split( '_' ) 
# 	if '_cru_' in fn:
# 		fn_parts = ['variable', 'metric', 'model_1', 'model_2', 'kind', 'month', 'year']
# 		fn_dict = dict( zip( fn_parts, fn_split ) )
# 		fn_dict.update( scenario='historical', model='cru_ts31' )
# 		# name_convention = [ 'variable', 'metric', 'model', 'scenario', 'experiment', 'begin_time', 'end_time' ]
# 	else:
# 		fn_parts = ['variable', 'metric', 'model', 'scenario', 'ensemble', 'month', 'year']
# 		fn_dict = dict( zip( fn_parts, fn_split ) )


# 	fn_switch = { 'cld':'_'.join([ 'cld','mean','pct','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif',
# 		'vap':'_'.join(['vap','mean','hPa','iem', fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif', 
# 		'tas':'_'.join(['tas','mean','C','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif',
# 		'hur':'_'.join(['hur','mean','pct','iem',fn_dict['model'],fn_dict['scenario'],fn_dict['month'], fn_dict['year'] ]) + '.tif' }

# 	output_filename = os.path.join( output_path, fn_switch[ fn_dict[ 'variable' ] ] )

# 	fn_list = fn.split( '.' )[0].split( '_' )
# 	return { i:j for i,j in zip( name_convention, fn_list )}
