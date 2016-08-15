# script to convert the newly generated Relative Humidity 
def convert_to_hur( tas_arr, vap_arr ):
	esa_arr = 6.112 * np.exp( 17.62 * tas_arr/ (243.12 + tas_arr) )
	# esa_arr = 6.112 * np.exp( 22.46 * tas_arr / (272.62 + tas_arr) )
	return vap_arr/esa_arr * 100
def convert_to_vap( tas_arr, hur_arr ):
	esa_arr = 6.112 * np.exp( 17.62 * tas_arr / (243.12 + tas_arr) )
	# esa_arr = 6.112 * np.exp( 22.46*tas_arr / (272.62 + tas_arr) )
	return (hur_arr * esa_arr) / 100
def run( x ):
	tas = rasterio.open( x[0] )
	hur = rasterio.open( x[1] )
	
	meta = tas.meta
	meta[ 'dtype' ] = 'float32' # set it to float32
	meta.update( compress='lzw' )
	meta.pop( 'transform' )

	tas_arr = tas.read( 1 )
	hur_arr = hur.read( 1 )
	vap_arr = convert_to_vap( tas_arr, hur_arr )

	# mask it:
	mask = tas.read_masks( 1 )
	vap_arr[ mask == 0 ] = tas.nodata

	# build an output filename from the input tas and write out
	output_filename = x[0].replace( 'tas', 'vap' )
	output_filename = output_filename.replace( '_C_', '_hPa_' )
	dirname = os.path.dirname( output_filename )

	if not os.path.exists( dirname ):
		os.makedirs( dirname )

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( vap_arr.astype( np.float32 ), 1 )
	return output_filename

if __name__ == '__main__':
	# import modules
	import os, glob, rasterio
	import numpy as np
	from pathos import multiprocessing as mp

	# read in the args
	tas_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323/tas/downscaled'
	hur_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323/hur/downscaled'
	tas_list = sorted( glob.glob( os.path.join( tas_path, 'tas*.tif' ) ) )
	hur_list = sorted( glob.glob( os.path.join( hur_path, 'hur*.tif' ) ) )

	# run in parallel
	pool = mp.Pool( 12 )
	out = pool.map( run, zip( tas_list, hur_list ) )
	pool.close()
