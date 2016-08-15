# script to bring the GIRR Radiation at top of Atmosphere through the clouds to generate 
def generate_nirr( cld_fn, girr_arr, output_filename ):
	'''
	generate nirr from girr and clouds for a common month

	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	** this is the c++ version given by Dave.  I am using R to calculate it here
	if (clds > -0.1) {
			nirr = cld.subset.v * (0.251 + (0.509*(1.0 - clds/100.0)))
		} else {
			nirr = -999.9; }
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	'''
	with np.errstate( over='ignore' ):
		cld = rasterio.open( cld_fn )
		mask = cld.read_masks( 1 )
		cld_arr = cld.read( 1 )
		nirr_arr = girr_arr * ( 0.251 + ( 0.509 * ( 1.0 - cld_arr/100 ) ) )
		nirr_arr[ mask == 0 ] = cld.nodata

	meta = cld.meta
	meta.update( compress='lzw' )
	if 'transform' in meta.keys():
		meta.pop( 'transform' )
	
	try:
		if not os.path.exists( os.path.dirname( output_filename ) ):
			os.makedirs( os.path.dirname( output_filename ) )
	except:
		pass

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( nirr_arr, 1 )
	return output_filename

if __name__ == '__main__':
	import os, glob, rasterio
	import numpy as np
	import pandas as pd
	from pathos.mp_map import mp_map

	# args
	ncores = 32

	# list the girr
	girr_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/girr_radiation_cmip3_process/IEM'
	girr = [ rasterio.open(i).read(1) for i in sorted( glob.glob( os.path.join( girr_dir, '*_fix.tif' ) ) ) ]
	months = ['01', '02', '03', '04', '05', '06', '07', '08', '09','10', '11', '12']
	# data dict for passing around in parallel
	girr = dict(zip( months,  girr ))

	# # list the cloud files for a series
	# input_path = '/Data/malindgren/cru_november_final/IEM/clouds/ar5'
	# 	# THERE IS AN ISSUE WITH THE IPSL...
	# models = [ 'MRI-CGCM3' ] # [ 'IPSL-CM5A-LR', 'GISS-E2-R' ] #[ 'MRI-CGCM3', 'CCSM4', 'GFDL-CM3' ] 
	# variables = [ 'cld' ] # we are only using clouds here
	# experiments = ['rcp26'] # [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	# path_list = [ os.path.join( input_path, model, variable, 'downscaled', experiment, '*.tif' ) for model in models for variable in variables for experiment in experiments ]

	# for path in path_list:
	# 	print( 'running: %s ' % path )
	# 	# path='/Data/malindgren/cru_november_final/IEM/ar5/MRI-CGCM3/cld/downscaled/*rcp26*.tif'
	# 	cld = pd.Series( glob.glob( path ) )
	# 	print( 'file count: %d' % len( cld ) )

	# 	output_filenames = cld.apply( lambda x: x.replace('cld', 'rsds').replace( '_pct_', '_MJ-m2-d1_' ) ).tolist()
	# 	month_grouper = cld.apply( lambda x: os.path.basename( x ).split( '.' )[0].split( '_' )[-2] )
	# 	args_list = [ ( cld, girr[ month ], out ) for cld, out, month in zip( cld.tolist(), output_filenames, month_grouper ) ]

	# 	# run it in parallel
	# 	out = mp_map( lambda x: generate_nirr( *x ), args_list, nproc=ncores )


	# # # CRU TS 3.23 Historical DATA VERSION
	# list the cloud files for a series
	path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/IEM/cru_ts31/cld/downscaled/*.tif'
	print( 'running: %s ' % path )
	cld = pd.Series( glob.glob( path ) )
	print( 'file count: %d' % len( cld ) )

	output_filenames = cld.apply( lambda x: x.replace('cld', 'rsds').replace( '_pct_', '_MJ-m2-d1_' ) ).tolist()
	month_grouper = cld.apply( lambda x: os.path.basename( x ).split( '.' )[0].split( '_' )[-2] )
	args_list = [ ( cld, girr[ month ], out ) for cld, out, month in zip( cld.tolist(), output_filenames, month_grouper ) ]

	# run it in parallel
	out = mp_map( lambda x: generate_nirr( *x ), args_list, nproc=ncores )

