# - - - - - - - - - - - - - - - - - - 
# masking NEW ALFRESCO Input Dataset
# - - - - - - - - - - - - - - - - - - 
def make_alfresco_compatible( fn, mask_fn, input_base_dir, output_base_dir ):
	import shutil, os, subprocess, rasterio
	import numpy as np	

	# split the input filename to parts
	dirname, basename = os.path.split( fn )
	
	# output filename to a new directory -- preserving the model/scenario/variable hierarchy
	out_fn = os.path.join( dirname.replace( input_base_dir, output_base_dir ), basename.replace( '_ar5_', '_alf_ar5_' ) )
	
	# read in the mask file where 0=nodata
	mask = rasterio.open( mask_fn ).read( 1 )
	
	try:
		# make sure the new dir exists...
		newdir = os.path.dirname( out_fn ) 
		if not os.path.exists( newdir ):
			os.makedirs( newdir )
	except:
		pass

	# copy mask_fn to new output filename
	_ = shutil.copy( mask_fn, out_fn )

	# gdalwarp -- subprocess... this is the best, easiest, fastest way... welcome to python.
	_ = subprocess.call([ 'gdalwarp', '-q', '-multi', '-dstnodata', 'None' ,'-srcnodata', 'None', fn, out_fn ])

	# read the array back in and mask it
	with rasterio.open( out_fn ) as rst:
		arr = rst.read( 1 )
		meta = rst.meta
	
	del rst

	# overwrite with a mask
	meta.update( nodata=-9999, crs={'init':'epsg:3338'}, dtype='float32', compress='lzw' )
	with rasterio.open( out_fn, 'w', **meta ) as out:
		arr[ mask == 0 ] = -9999
		out.write( arr.astype( np.float32 ), 1 )

	return out_fn

if __name__ == '__main__':
	import rasterio, subprocess, os
	import numpy as np
	# from pathos.mp_map import mp_map
	from functools import partial
	import multiprocessing as mp
	import argparse

	# # # parse the commandline arguments
	# parser = argparse.ArgumentParser( description='reformat AR5-CMIP5 data to 1km resolution for ALFRESCO' )
	# parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="cmip5 model name (exact)" )
	# parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="cmip5 scenario name (exact)" )
	# parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	# args = parser.parse_args()

	# # unpack args
	# model = args.model
	# scenario = args.scenario
	# variable = args.variable

	# # # # STEP1 MAKE A MASK FROM THE Version 1 CMIP3 DATA
	input_base_dir = os.path.join( '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled', model, scenario, variable )
	output_base_dir = os.path.join( '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/alfresco_1km', model, scenario, variable )

	# HARDWIRED ARGS...
	alf_v1 = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/alfresco_formatting/tas_mean_C_alf_ar4_cccma_cgcm3_1_sresa2_01_2001.tif'
	mask_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/alfresco_formatting/alf_v1_mask_1km.tif'

	# make the mask file --- this only need be run once...
	# # this is the old data
	# rst = rasterio.open( alf_v1 )
	# mask = rst.read_masks(1)
	# mask[ mask == 255 ] = 1
	# meta = rst.meta
	# del rst

	# # write the extracted mask to disk.
	# meta.update( crs={'init':'epsg:3338'}, compress='lzw', dtype='float32' )
	# with rasterio.open( mask_fn, 'w', **meta ) as out:
	# 	out.write( mask.astype( np.float32 ), 1 )

	# del out

	# STEP2 list all the data and make args for passing to the parallel function
	filelist = [ os.path.join( r, fn ) for r,s,files in os.walk( input_base_dir ) for fn in files if fn.endswith( '.tif' ) ]
	f = partial( make_alfresco_compatible, mask_fn=mask_fn, input_base_dir=input_base_dir, output_base_dir=output_base_dir )

	# # # # STEP3 run in parallel
	pool = mp.Pool( 32 )
	done = pool.map( f, filelist )
	# done = mp_map( f, filelist, nproc=32 )
	pool.close()
	pool.join()

	pool = None
	del pool
