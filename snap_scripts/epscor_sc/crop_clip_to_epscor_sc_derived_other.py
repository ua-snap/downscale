def make_args( rst_fn, shp_fn, output_path ):
	'''
	generate tuples of the arguments needed to run crop_clip with all data
	'''
	import os
	dirname, filename = os.path.split( rst_fn )
	# ['dof', '5ModelAvg', 'rcp85', '2090', '2099']
	variable, model, scenario, begin_year, end_year = filename.split('.')[0].split('_')
	# variable, metric, units, project, model, scenario, month, year = filename.split('.')[0].split('_')
	# overcome a naming error from someone else
	if model == '5modelAvg':
		model = '5ModelAvg'
		filename = filename.replace( '5modelAvg', '5ModelAvg' )
	out_fn = os.path.join( output_path, model, scenario, variable, filename )
	out_fn = out_fn.replace( '_c_', '_C_' )
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

	# proj4string = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
	proj4string = 'EPSG:3338'
	subprocess.call(['gdalwarp', '-q', '-t_srs', proj4string, '-tap','-overwrite', \
						'-co', 'COMPRESS=LZW','-tr', '2000', '2000', \
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
	from pathos.mp_map import mp_map

	# setup args
	base_path = '/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models'
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/EPSCOR_SC_DELIVERY_AUG2016/derived/grids/monthly_decadals'
	ncpus = 32
	subdomain_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/Kenai_StudyArea.shp'
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4', '5ModelAvg' ]

	# list up all the args we want to run through the multicore clipping
	fn_list = []
	for root, subs, files in os.walk( base_path ):
		if any( [ model in root for model in models ] ) == True:
			if 'derived' in root:
				if len( [ fn for fn in files if fn.endswith( '.tif' ) ] ) > 0:
					fn_list = fn_list + glob.glob( os.path.join( root, '*.tif' ) )

	args_list = [ make_args( rst_fn, subdomain_fn, output_path ) for rst_fn in fn_list if 'dof' in rst_fn or 'dot' in rst_fn or 'logs' in rst_fn ]
	out = mp_map( lambda x: wrap( x ), args_list, nproc=32 )
	