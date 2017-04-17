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

	# proj4string = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
	proj4string = 'EPSG:3338'
	subprocess.call(['gdalwarp', '-q', '-tap','-overwrite', '-tap' ,'-t_srs', proj4string,'-co', 'COMPRESS=LZW', '-tr', '2000', '2000', 
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
	from pathos.mp_map import mp_map

	# setup args
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/decadal_monthlies'
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/EPSCOR_SC_DELIVERY_NOV2016/derived/grids/decadal_monthlies'
	ncpus = 32
	subdomain_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/Kenai_StudyArea.shp'
	# list up all the args we want to run through the multicore clipping
	args_list = []
	for root, subs, files in os.walk( base_path ):
		if len( [ fn for fn in files if fn.endswith( '.tif' ) if 'dof_mean_decadal' in fn or 'dot_mean_decadal' in fn or 'logs_mean_decadal' in fn ] ) > 0:
			args_list = args_list + [ ( subdomain_fn, os.path.join( root, fn ), os.path.join( root, fn ).replace( base_path, output_path ) ) for fn in files ]
	
	out = mp_map( wrap, args_list, nproc=ncpus )
