# script to convert the newly generated Relative Humidity 
def convert_to_hur( tas_arr, vap_arr ):
	esa_arr = 6.112 * np.exp( 17.62 * tas_arr/ (243.12 + tas_arr) )
	# esa_arr = 6.112 * np.exp( 22.46 * tas_arr / (272.62 + tas_arr) )
	return vap_arr/esa_arr * 100
def convert_to_vap( tas_arr, hur_arr ):
	''' create relative humidity from the CRU tas / vap '''
	esa_arr = 6.112 * np.exp( 17.62 * tas_arr/ (243.12 + tas_arr) )
	# esa_arr = 6.112 * np.exp( 22.46 * tas_arr / (272.62 + tas_arr) )
	return (hur_arr*esa_arr)/100
def run( x ):
	tas = rasterio.open( x[0] )
	hur = rasterio.open( x[1] )
	mask = tas.read_masks( 1 ) != 0
	
	meta = tas.meta
	meta[ 'dtype' ] = 'float32' # set it to float32
	meta.update( compress='lzw' )
	meta.pop( 'transform' )

	tas_arr = tas.read( 1 )
	hur_arr = hur.read( 1 )
	vap_arr = np.copy( tas_arr )
	vap_arr[ mask ] = convert_to_hur( tas_arr[ mask], hur_arr[ mask ] )

	# filename-fu
	output_filename = x[0].replace( 'tas', 'vap' ).replace( '_C_', '_hPa_' )
	dirname = os.path.dirname( output_filename )

	if not os.path.exists( dirname ):
		os.makedirs( dirname )

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( np.around( vap_arr.astype( np.float32 ), 2 ), 1 )
	return output_filename

if __name__ == '__main__':
	import os, glob, rasterio
	import numpy as np
	from pathos.mp_map import mp_map

	# read in the args
	base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016'
	for model in ['CCSM4', 'GFDL-CM3', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R']:
		for scenario in ['historical', 'rcp26','rcp45','rcp60','rcp85']:
			tas_list = sorted( glob.glob( os.path.join( base_path, 'downscaled', model, scenario, 'tas', '*.tif' ) ) )
			hur_list = sorted( glob.glob( os.path.join( base_path, 'downscaled', model, scenario, 'hur', '*.tif' ) ) )

			# make args to pass to the run function
			args = zip( tas_list, hur_list )
			# run in parallel
			out = mp_map( run, args, nproc=32 )
