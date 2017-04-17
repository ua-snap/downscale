# FORCE ROUNDING / TRUNCATION
# ---------------------------

def run( fn ):
	import rasterio
	import numpy as np
	import os

	dirname, basename = os.path.split( fn )
	variable = basename.split( '_' )[0] # get var from the fn first elem

	# open in update mode
	with rasterio.open( fn ) as rst:
		arr = rst.read( 1 )
		meta = rst.meta
		meta.update( compress='lzw', nodata=-9999 )
		
		try:
			meta.pop( 'transform' )
		except:
			pass

		if variable == 'pr':
			arr[ arr != -9999 ] = np.rint( arr[ arr != -9999 ] )
		else:
			arr[ arr != -9999 ] = np.round( arr[ arr != -9999 ], 1 )

	with rasterio.open( fn, 'w', **meta ) as out:
		out.write( arr, 1 )

if __name__ == '__main__':
	import os
	from pathos.mp_map import mp_map

	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled'
	files = [ os.path.join( root, fn ) for root, subs, files in os.walk( base_path ) for fn in files if fn.endswith( '.tif' ) and 'tas' in fn or 'pr_' in fn or 'swi_' in fn ]

	done = mp_map( run, files, nproc=32 )

