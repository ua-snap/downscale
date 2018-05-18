# convert sunp to clt
import os, glob, shutil, rasterio
import numpy as np

base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/insolation_L48/climatologies'
variable = 'sunp'
out_variable = 'clt'

# list data
files = glob.glob( os.path.join( base_path, variable, '*'+variable+'*.tif' ) )

for fn in files:
	print( fn )
	dirname, basename = os.path.split( fn )
	new_dir = dirname.replace( variable, out_variable )
	if not os.path.exists( new_dir ):
		os.makedirs( new_dir )
	
	out_fn = fn.replace( variable, out_variable )

	# open the raster
	with rasterio.open( fn ) as rst:
		sunp_arr = rst.read( 1, masked=True )
		sunp_arr.fill_value = rst.nodata
		meta = rst.meta
		meta.update( compress='lzw' )	

	clt_arr = np.around( (100.0 - sunp_arr).astype( np.float32 ), 1 )

	with rasterio.open( out_fn, 'w', **meta ) as out:
		out.write( clt_arr.filled(), 1 )

