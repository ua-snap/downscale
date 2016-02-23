# script to convert the newly generated Relative Humidity 
def convert_to_hur( tas_arr, vap_arr ):
	import numpy as np
	with np.errstate( over='ignore' ):
		esa_arr = 6.112 * np.exp( 17.62 * tas_arr/ (243.12 + tas_arr) )
		# esa_arr = 6.112 * np.exp( 22.46 * tas_arr / (272.62 + tas_arr) )
		return vap_arr/esa_arr * 100
def convert_to_vap( tas_arr, hur_arr ):
	import numpy as np
	with np.errstate( over='ignore' ):
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

	# build an output filename from the input tas and write out -- changed to deal with pathing!
	output_filename = x[1].replace( 'hur', 'vap' )
	output_filename = output_filename.replace( '_metric_', '_hPa_' )
	# output_filename = x[0].replace( 'tas', 'vap' )
	# output_filename = output_filename.replace( '_C_', '_hPa_' )
	dirname = os.path.dirname( output_filename )
	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( vap_arr.astype( np.float32 ), 1 )
	return output_filename

if __name__ == '__main__':
	# import modules
	import os, glob, rasterio
	import numpy as np
	from pathos import multiprocessing as mp

	# args
	ncores = 40
	tas_input_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_november_final/ar5'
	hur_input_path = '/Data/malindgren/cru_november_final/ar5'

	models = [ 'IPSL-CM5A-LR', 'GISS-E2-R', 'MRI-CGCM3', 'CCSM4', 'GFDL-CM3' ]

	for model in models:
		print model
		tas_files = sorted( glob.glob( os.path.join( tas_input_path, model, 'tas', 'downscaled', '*.tif' ) ) )
		hur_files = sorted( glob.glob( os.path.join( hur_input_path, model, 'hur', 'downscaled', '*.tif' ) ) )

		# combine the sorted lists which should now be in a common order...
		tas_hur_list = zip( tas_files, hur_files )

		# run in parallel
		pool = mp.Pool( processes=ncores )
		out = pool.map( run, tas_hur_list )
		pool.close()



# def return_files( input_path, var ):
# 	output_files = []
# 	for root, subs, files in os.walk( input_path ):
# 		# print root
# 		if root.endswith( 'downscaled' ) and len( files ) != 0 and var in root:
# 			pool = mp.Pool( processes=ncores )
# 			files = pool.map( lambda x: os.path.join( root, x ), files )
# 			pool.close()
# 			output_files.append( files )
# 	return output_files

