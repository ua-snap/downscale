def make_args( rst_fn, shp_fn, output_path ):
	'''
	generate tuples of the arguments needed to run crop_clip with all data
	'''
	import os
	dirname, filename = os.path.split( rst_fn )
	variable, metric, units, project, model, scenario, month, year = filename.split('.')[0].split('_')
	# overcome a naming error from someone else
	if model == '5modelAvg':
		model = '5ModelAvg'
		filename = filename.replace( '5modelAvg', '5ModelAvg' )
	out_fn = os.path.join( output_path, model, scenario, variable, filename )
	return ( shp_fn, rst_fn, out_fn )

def crop_clip( shp_fn, rst_fn, out_fn ):
	'''
	crop/clip to the shapefile we want using gdalwarp.
	'''
	import subprocess, os
	try:
		dirname = os.path.dirname( out_fn )
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass

	subprocess.call(['gdalwarp', '-q', '-tap','-overwrite', \
						'-co', 'COMPRESS=LZW', '-tr', '2000', '2000', \
						'-dstnodata', '-3.4e+38','-crop_to_cutline', '-cutline', 
						shp_fn, rst_fn, out_fn ])
	return out_fn

def wrap( x ):
	''' wrapper for clean multiprocessing call to pool.map '''
	return crop_clip( *x )

if __name__ == '__main__':
	import os, glob, itertools, rasterio
	import xarray as xr
	import numpy as np
	import pandas as pd
	from pathos import multiprocessing as mp

	# setup args
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_v2/CRU_TS323'
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_v2_clipped'
	ncpus = 32
	subdomain_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/Kenai_StudyArea.shp'

	# list up all the args we want to run through the multicore clipping
	fn_list = []
	for root, subs, files in os.walk( base_path ):
		if len( [ fn for fn in files if fn.endswith( '.tif' ) ] ) > 0:
			fn_list = fn_list + glob.glob( os.path.join( root, '*.tif' ) )

	args_list = [ make_args( rst_fn, subdomain_fn, output_path ) for rst_fn in fn_list ]

	pool = mp.Pool( ncpus )
	out = pool.map( wrap, args_list )
	pool.close()
	pool.join()
	pool.terminate()
