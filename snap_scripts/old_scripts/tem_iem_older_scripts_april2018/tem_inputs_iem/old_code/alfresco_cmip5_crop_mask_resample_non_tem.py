def resample_to_1km( x ):
	'''
	template_raster_mask should be a mask in in the res/extent/origin/crs of the 
	existing TEM IEM products.
	'''
	fn, output_filename, template_raster_mask_fn = x
	import rasterio, os
	from rasterio.warp import RESAMPLING, reproject
	import numpy as np

	rst = rasterio.open( fn )
	rst_arr = rst.read( 1 )

	fn = os.path.basename( fn )
	fn_split = fn.split( '.' )[0].split( '_' ) 

	template_raster_mask = rasterio.open( template_raster_mask_fn )
	template_arr = template_raster_mask.read( 1 )
	template_meta = template_raster_mask.meta
	template_meta.update( compress='lzw', nodata=rst.nodata, dtype='float32' )
	
	if 'transform' in template_meta.keys():
		template_meta.pop( 'transform' )

	output_arr = np.empty_like( template_arr.astype( np.float32 ) )
	output_arr[ template_arr == 0 ] = rst.nodata

	epsg3338 = {'init':'epsg:3338'}
	reproject( rst_arr, output_arr, src_transform=rst.affine, src_crs=epsg3338, src_nodata=rst.nodata, \
			dst_transform=template_raster_mask.affine, dst_crs=epsg3338,\
			dst_nodata=rst.nodata, resampling=RESAMPLING.cubic_spline ) # , num_threads=2

	# REQUIRES UPDATING! 
	# this is where we would need to truncate the data to 2 decimal places for the tas
	#  and we need have a dtype diff I think for pr = int32 and tas = float32
	# this needs to match the AR4 data! check the CKAN_Data data.

	# get indexes of non-masked values
	ind = np.where( output_arr >= 0 )
	if 'pr_' in fn:
		# truncate the values to 0 decimal places
		new_vals = np.ceil( output_arr[ ind ] )
	elif 'tas_' in fn:
		# round to 2 decimal places
		new_vals = np.round( output_arr[ ind ], 2 )
	else:
		AttributeError( 'tool only built to work with tas / pr at this time.' )

	# pass the prepped values to the output_arr
	output_arr[ ind ] = new_vals

	with rasterio.open( output_filename, 'w', **template_meta ) as out:
		output_arr[ template_arr == 0 ] = rst.nodata
		out.write( output_arr, 1 )
	return output_filename

def path_slicer_ar5( input_path ):
	'''
	slicer works for the /Data AR5 data outputs produced by @leonawicz
	and is hardwired to that data structure.
	'''
	scenario, model, variable = input_path.split( os.path.sep )[-3:]
	return scenario, model, variable

def generate_output_fn( input_fn, output_path, name ):
	'''
	add id as an element in the new filename
	'''
	base = os.path.basename( input_fn )
	basesplit = base.split( '_' )
	basesplit.insert( 3, name ) #name = 'iem' or 'akcan'
	return os.path.join( output_path, '_'.join( basesplit ) )

if __name__ == '__main__':
	import os, glob, rasterio
	import numpy as np1
	import pandas as pd
	from functools import partial
	import itertools
	from pathos import multiprocessing as mp
	from pathos.mp_map import mp_map
	import argparse

	parser = argparse.ArgumentParser( description='crop mask and resample the alf inputs to different extents and resolutions' )
	parser.add_argument( '-p', '--input_path', action='store', dest='input_path', type=str, help='path to input files' )
	parser.add_argument( '-o', '--output_base_path', action='store', dest='output_base_path', type=str, help='path to output base directory' )
	parser.add_argument( '-m', '--template_raster_mask_fn', action='store', dest='template_raster_mask_fn', type=str, help='path to mask raster' )
	parser.add_argument( '-g', '--group', action='store', dest='group', type=str, help='name to identify raster type in filename. i.e. "iem", "alf" ')

	args = parser.parse_args()
	input_path = args.input_path
	output_base_path = args.output_base_path
	template_raster_mask_fn = args.template_raster_mask_fn
	group = args.group
	ncores = 32

	# more prep
	scenario, model, variable = path_slicer_ar5( input_path )
	template_raster_mask = rasterio.open( template_raster_mask_fn )	
	output_path = os.path.join( output_base_path, scenario, model, variable )

	try:
		if not os.path.exists( output_path ):
			os.makedirs( output_path )
	except:
		pass

	file_list = glob.glob( os.path.join( input_path, '*.tif' ) )
	output_filenames = [ generate_output_fn( input_fn, output_path, group ) for input_fn in file_list ]

	args = zip( file_list, output_filenames, itertools.repeat( template_raster_mask_fn, len(file_list) ) ) 
	_ = mp_map( resample_to_1km, args, nproc=ncores )


# # # # RUN THE ABOVE # # # #
# import os, glob
# import numpy as np

# # list the data we want
# input_path = '/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models'
# out = [ root for root, subs, files in os.walk( input_path ) \
# 		if len( glob.glob( os.path.join( root, '*.tif' ) ) ) > 0 and not 'derived' in root ]
# input_paths = np.unique( out ).tolist()

# models = [ 'IPSL-CM5A-LR', 'MPI-ESM-LR', 'MRI-CGCM3', '5modelAvg' ] # 'CCSM4', 'CNRM-CM5', 'GFDL-CM3', 'GISS-E2-R']
# scenarios = [ 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]

# # # # IEM 1km
# # output_base_path = '/atlas_scratch/malindgren/CMIP5/IEM_1km'
# # template_raster_mask_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/extents/IEM_Mask_1km.tif'
# # group = 'iem'

# # for input_path in input_paths:
# # 	if model in input_path:
# # 		print input_path
# # 		os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# # 		os.system( 'python alfresco_cmip5_crop_mask_resample_non_tem.py -p ' + input_path + ' -o ' + output_base_path + ' -m ' + template_raster_mask_fn + ' -g ' + group )

# # # AKCAN 1km
# output_base_path = '/atlas_scratch/malindgren/CMIP5/AKCAN_1km'
# template_raster_mask_fn = '/workspace/Shared/Users/malindgren/test_stuff/new_mask.tif'
# group = 'alf'

# for input_path in input_paths:
# 	for model in models:
# 		if model in input_path:
# 			print input_path
# 			os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# 			os.system( 'python alfresco_cmip5_crop_mask_resample_non_tem.py -p ' + input_path + ' -o ' + output_base_path + ' -m ' + template_raster_mask_fn + ' -g ' + group )

# # # CRU DATA RUN EXAMPLE:
# # # RUN THE ABOVE # # # #
# import os, glob
# import numpy as np

# # list the data we want
# input_path = '/Data/Base_Data/Climate/AK_CAN_2km/historical/CRU/CRU_TS32'
# out = [ root for root, subs, files in os.walk( input_path ) \
# 		if len( glob.glob( os.path.join( root, '*.tif' ) ) ) > 0 and not 'derived' in root ]
# input_paths = np.unique( out ).tolist()

# # # # IEM 1km
# # output_base_path = '/atlas_scratch/malindgren/CMIP5/IEM_1km'
# # template_raster_mask_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/extents/IEM_Mask_1km.tif'
# # group = 'iem'

# # for input_path in input_paths:
# # 	if model in input_path:
# # 		print input_path
# # 		os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# # 		os.system( 'python alfresco_cmip5_crop_mask_resample_non_tem.py -p ' + input_path + ' -o ' + output_base_path + ' -m ' + template_raster_mask_fn + ' -g ' + group )

# # # AKCAN 1km
# output_base_path = '/atlas_scratch/malindgren/CMIP5/AKCAN_1km'
# template_raster_mask_fn = '/workspace/Shared/Users/malindgren/test_stuff/new_mask.tif'
# group = 'alf'

# for input_path in input_paths:
# 	print input_path
# 	os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# 	os.system( 'python alfresco_cmip5_crop_mask_resample_non_tem.py -p ' + input_path + ' -o ' + output_base_path + \
# 					' -m ' + template_raster_mask_fn + ' -g ' + group )



# # # NOTES: 
# some setup:
# cru_path = '/Data/Base_Data/Climate/AK_CAN_2km/historical/CRU/CRU_TS32'
# ar5_path = '/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models'
# input_path = '/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models'
# output_base_path = '/atlas_scratch/malindgren/CMIP5'
