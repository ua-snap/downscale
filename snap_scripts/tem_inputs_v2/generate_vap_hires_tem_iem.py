# # # # # # # # # # # # # # 
# # convert tas/hur to vap
# # # # # # # # # # # # # # 

def convert_to_vap( tas_arr, hur_arr ):
	''' create relative humidity from the CRU tas / vap '''
	esa_arr = 6.112 * np.exp( 17.62 * tas_arr/ (243.12 + tas_arr) )
	# esa_arr = 6.112 * np.exp( 22.46 * tas_arr / (272.62 + tas_arr) )
	return (hur_arr*esa_arr)/100

def make_vap( hur_fn, tas_fn ):
	''' make vapor pressure from hur and tas.'''
	hur = rasterio.open( hur_fn )
	tas = rasterio.open( tas_fn )

	hur_arr = hur.read( 1 )
	tas_arr = tas.read( 1 )
	with np.errstate( all='ignore' ):
		mask = hur.read_masks( 1 )
		vap_arr = convert_to_vap( tas_arr, hur_arr )
		
	vap_arr[ mask == 0 ] = hur_arr.min()
	vap_arr = np.around( vap_arr, 1 ) # roundit
	vap_arr[ mask == 0 ] = hur_arr.min() # reset mask

	output_filename = hur_fn.replace( 'hur', 'vap' ).replace( 'pct', 'hPa' )
	dirname = os.path.dirname( hur_fn ).replace( 'hur', 'vap' )

	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass 

	meta = hur.meta
	meta.update( compress='lzw' )
	_ = meta.pop( 'transform' )

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( vap_arr, 1 )

	return output_filename

def wrap_make_vap( x ):
	''' runner for multiprocessing '''
	return make_vap( *x )

if __name__ == '__main__':
	import os, rasterio, itertools, functools, glob
	import numpy as np
	from pathos.mp_map import mp_map

	# # vars
	hur_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/downscaled'
	tas_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled'
	# models =  [ 'IPSL-CM5A-LR', 'GISS-E2-R', 'MRI-CGCM3', 'NCAR-CCSM4', 'GFDL-CM3' ]
	# scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	# # cru run
	models = [ 'CRU_TS323' ]
	scenarios = [ 'historical' ]

	out = []
	for model, scenario in itertools.product( models, scenarios ):
		# list files
		hur_files = glob.glob( os.path.join( hur_dir, model, scenario, 'hur', '*.tif' ) )
		basenames = [ os.path.basename( fn ).replace('hur_mean_pct', 'tas_mean_C') for fn in hur_files ]
		tas_files = [ os.path.join( tas_dir, model, scenario, 'tas', fn.replace( 'cru_ts323', 'CRU_TS323') ) for fn in basenames ]
		files = zip( hur_files, tas_files )
		
		# run it in parallel
		out = out + [ mp_map( wrap_make_vap, files, nproc=32 ) ]
