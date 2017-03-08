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
	subprocess.call(['gdalwarp', '-r', 'bilinear', '-q', '-overwrite','-t_srs', proj4string, '-co', 'COMPRESS=LZW', 
						'-tr', '1000', '1000','-srcnodata', '-9999', '-dstnodata', '-9999', 
						'-cutline', shp_fn, '-crop_to_cutline', rst_fn, out_fn ])
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
	base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/downscaled_1km'
	output_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/downscaled_alfresco'
	ncpus = 32
	# subdomain_fn = '/workspace/Shared/Tech_Projects/Alaska_IEM/project_data/ALFRESCO/aiem_domain_fixed/AIEM_domain.shp'
	subdomain_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/masks/alfresco_akcan_mask.shp'
 
	# list up all the args we want to run through the multicore clipping
	args_list = []
	for root, subs, files in os.walk( base_path ):
		tif_files = [ fn for fn in files if fn.endswith( '.tif' ) and 'tas_' in fn ]
		# tif_files = [ fn for fn in files if fn.endswith( '.tif' ) and 'pr_' in fn ]
		# tif_files = [ fn for fn in files if fn.endswith( '.tif' ) and 'vap_' in fn ]
		# tif_files = [ fn for fn in files if fn.endswith( '.tif' ) and 'rsds_' in fn ]
		
		if len( tif_files ) > 0:
			args_list = args_list + [ ( subdomain_fn, os.path.join( root, fn ), os.path.join( root, fn ).replace( base_path, output_path ) ) for fn in tif_files ]
	
	out = mp_map( wrap, args_list, nproc=ncpus )
