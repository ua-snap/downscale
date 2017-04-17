# standardize the SNAP data to the new standard
# oob = -9999 crs = 3338 change cru namings...
def list_all_data( base_path ):
	return [os.path.join(root, fn) for root, subs, files in os.walk( base_path ) for fn in files if fn.endswith( '.tif' ) ]

def fix_naming( fn ):
	dirname, basename = os.path.split( fn )
	basename = basename.replace( 'cru_ts323', 'CRU_TS323' ).replace( 'cru_TS323', 'CRU_TS323' )
	basename = basename.replace( 'cru_ts40', 'CRU_TS40' ).replace( 'cru_TS40', 'CRU_TS40' )
	basename = basename.replace( '5MODELAVG', '5ModelAvg' )
	return os.path.join( dirname, basename )

def update_crs_oob( fn, base_path, output_path ):
	# reproject to 3338, update nodata value and make sure its compressed
	# I have chosen not to do this since it is off by a precision we cannot even be certain of given our
	# data constraints.  I am just going to set it to 3338 without a reproject
	# os.system( 'gdalwarp -overwrite -t_srs EPSG:3338 -tr 2000 2000 -dstnodata -9999 -co COMPRESS=LZW {} {}'.format( fn, out_fn ) )
	with rasterio.open( fn ) as rst:
		arr = rst.read( 1 )
		meta = rst.meta
		meta.update( compress='lzw', count=1, nodata=-9999, crs={'init':'EPSG:3338'} )
		old_nodata = rst.nodata
		
	if 'transform' in meta.keys():
		_ = meta.pop( 'transform' )


	# overwrite the existing file...
	with rasterio.open( fn, 'w', **meta ) as out:
		arr[ arr == -3.39999995e+38 ] = -9999.0
		out.write( arr, 1 )

	# rename the existing file to the proper naming convention
	output_filename = fix_naming( fn.replace( base_path, output_path ) )

	# make sure the new directory actually exists
	dirname, basename = os.path.split( output_filename )
	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass
	
	# rename the newly generated and closed file
	os.rename( fn, output_filename )
	return output_filename

if __name__ == '__main__':
	import numpy as np
	import rasterio, os
	from pathos.mp_map import mp_map
	from functools import partial

	# base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids' # downscaled' 
	# output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids' # downscaled'

	# base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/downscaled'
	# output_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/downscaled' 

	# can be the same directory with an overwrite...
	base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'

	# list the data
	files = list_all_data( base_path )
	f = partial( update_crs_oob, base_path=base_path, output_path=output_path )
	done = mp_map( f, files, nproc=32 )
