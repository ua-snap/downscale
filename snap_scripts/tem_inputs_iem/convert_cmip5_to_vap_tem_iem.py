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
	mask = tas.read_masks( 1 ) != 0
	
	meta = tas.meta
	meta[ 'dtype' ] = 'float32' # set it to float32
	meta.update( compress='lzw' )
	meta.pop( 'transform' )

	tas_arr = tas.read( 1 )
	hur_arr = hur.read( 1 )
	vap_arr = np.copy( tas_arr )
	vap_arr[ mask ] = convert_to_vap( tas_arr[ mask], hur_arr[ mask ] )

	# filename-fu
	output_filename = x[0].replace( 'tas', 'vap' ).replace( '_C_', '_hPa_' )
	dirname = os.path.dirname( output_filename )

	if not os.path.exists( dirname ):
		os.makedirs( dirname )

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( vap_arr.astype( np.float32 ), 1 )
	return output_filename

if __name__ == '__main__':
	import os, glob, rasterio
	import numpy as np
	import xarray as xr
	from pathos.mp_map import mp_map

	# read in the args
	base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016'
	cru_path = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/'
	args = []
	for model in ['CCSM4', 'GFDL-CM3', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R']:
		for scenario in ['historical', 'rcp26','rcp45','rcp60','rcp85']:
			tas_list = sorted( glob.glob( os.path.join( base_path, 'downscaled', model, scenario, 'tas', '*.tif' ) ) )
			hur_list = sorted( glob.glob( os.path.join( base_path, 'downscaled', model, scenario, 'hur', '*.tif' ) ) )

			# make args to pass to the run function
			args = args + zip( tas_list, hur_list )
	
	# run in parallel
	out = mp_map( run, args, nproc=32 )

	# # # CONVERT CL20 2km to vap
	tas_list = sorted(glob.glob( os.path.join( base_path, 'cru', 'cru_cl20', 'tas', '*.tif' ) ))
	hur_list = sorted(glob.glob( os.path.join( base_path, 'cru', 'cru_cl20', 'hur', '*.tif' ) ))
	args = zip( tas_list, hur_list )
	out = mp_map( run, args, nproc=12 )

	# # # CONVERT CRU TS323 vap/tas to hur --> output to a non CF-compliant NetCDF that will be read back in with xarray
	tas = xr.open_dataset( '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.tmp.dat.nc' )
	vap = xr.open_dataset( '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.vap.dat.nc' )
	hur = convert_to_hur( tas.tmp, vap.vap )
	hur_ds = hur.to_dataset( name='hur' )
	hur_ds.to_netcdf( '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.hur.SNAP_derived.dat.nc' )

