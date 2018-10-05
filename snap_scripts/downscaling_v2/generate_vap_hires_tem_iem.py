# # # # # # # # # # # # # # 
# # convert tas/hurs to vap
# # # # # # # # # # # # # # 

def convert_to_vap( tas_arr, hur_arr ):
	''' create relative humidity from the CRU tas / vap '''
	esa_arr = 6.112 * np.exp( 17.62 * tas_arr/ (243.12 + tas_arr) )
	# esa_arr = 6.112 * np.exp( 22.46 * tas_arr / (272.62 + tas_arr) )
	return (hur_arr*esa_arr)/100

def make_vap( hur_fn, tas_fn, out_fn ):
	''' make vapor pressure from hur and tas.'''
	with rasterio.open( hur_fn ) as hur:
		hur_arr = hur.read( 1 )
		mask = hur.read_masks( 1 )
		nodata = hur.nodata
		meta = hur.meta
		meta.update( compress='lzw' )
	
	with rasterio.open( tas_fn ) as tas:
		tas_arr = tas.read( 1 )

	with np.errstate( all='ignore' ):
		vap_arr = convert_to_vap( tas_arr, hur_arr )
		
	vap_arr = np.around( vap_arr, 2 ) # roundit
	vap_arr[ mask == 0 ] = nodata # reset mask

	dirname, basename = os.path.split( out_fn )
	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass 

	with rasterio.open( out_fn, 'w', **meta ) as out:
		out.write( vap_arr, 1 )

	return out_fn

def wrap_make_vap( x ):
	''' runner for multiprocessing '''
	return make_vap( *x )

if __name__ == '__main__':
	import os, rasterio, itertools, functools, glob
	import numpy as np
	import multiprocessing as mp
	
	# # vars
	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	ncpus = 64
	
	# list ALL relative humidity
	hurs_files = [ os.path.join(r,fn) for r,s,files in os.walk( base_dir ) for fn in files if fn.endswith( '.tif' ) and 'hurs_' in fn and '_anom.tif' not in fn and 'CRU' in fn ]

	# since the pathing is the same except for variable, metric, units we can just change the list to make a tas list
	tas_files = [ fn.replace('/hurs','/tas').replace('mean_pct','mean_C') for fn in hurs_files ]

	# make the output files from one of the lists
	output_filenames = [ fn.replace('/hurs','/vap').replace('mean_pct','mean_hPa') for fn in hurs_files ]

	# run processing in parallel
	args = zip( hurs_files, tas_files, output_filenames )
	pool = mp.Pool( ncpus )
	out = pool.map( wrap_make_vap, args )
	pool.close()
	pool.join()