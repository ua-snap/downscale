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

	subprocess.call(['gdalwarp', '-q', '-tap','-overwrite', '-co', 'COMPRESS=LZW', '-tr', '2000', '2000', 
						'-srcnodata', '-3.4e+38', '-dstnodata', '-3.4e+38', '-crop_to_cutline', '-cutline', 
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
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_v2'
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_v2_clipped'
	ncpus = 32
	subdomain_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/Kenai_StudyArea.shp'

	# list up all the args we want to run through the multicore clipping
	args_list = []
	for root, subs, files in os.walk( base_path ):
		tif_files = [ fn for fn in files if fn.endswith( '.tif' ) ]
		if len( tif_files ) > 0:
			args_list = args_list + [ ( subdomain_fn, os.path.join( root, fn ), os.path.join( root, fn ).replace( base_path, output_path ) ) for fn in tif_files ]
	
	pool = mp.Pool( ncpus )
	out = pool.map( wrap, args_list )
	pool.close()
	pool.join()
	pool.terminate()