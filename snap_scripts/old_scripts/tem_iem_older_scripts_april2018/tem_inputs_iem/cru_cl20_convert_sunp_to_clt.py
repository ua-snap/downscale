# convert sunp to clt
import os, glob, shutil, rasterio
import numpy as np

base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/cru/cru_cl20'
variables = ['sunp']
new_names = ['clt']

months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
months_lookup = { count+1:month for count, month in enumerate( months ) }

for variable, new_name in zip(variables, new_names):
	files = glob.glob( os.path.join( base_path, variable, '*'+variable+'*.tif' ) )
	
	for fn in files:
		print( fn )
		dirname, basename = os.path.split( fn )
		new_dir = dirname.replace( variable, new_name )
		if not os.path.exists( new_dir ):
			os.makedirs( new_dir )
		
		# fix single digit months
		mon = basename.split('_')[-2]
		basename = basename.replace( '_'+mon+'_', '_'+months_lookup[ int(mon) ]+'_' )

		# open the raster
		rst = rasterio.open( fn )
		arr = rst.read( 1, masked=True )
		arr.fill_value = rst.nodata
		clt_arr = np.round( (100.0 - arr).astype( np.float32 ), 1 )

		meta = rst.meta
		meta.update( compress='lzw' )		

		new_fn = os.path.join( new_dir, basename.replace( variable, new_name ) )
		with rasterio.open( new_fn, 'w', **meta ) as out:
			out.write( clt_arr.filled(), 1 )


		# change variable naming
		# shutil.copy( fn, new_fn )
