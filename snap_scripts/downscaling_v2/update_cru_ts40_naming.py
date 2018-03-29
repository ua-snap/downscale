def rename_cru_ts40( fn ):
	''' simple rename of the cru_ts40 data to the more proper naming convention. '''
	out_fn = fn.replace( 'cru_ts40', 'CRU-TS40' )
	return os.rename( fn, out_fn )
	
if __name__ == '__main__':
	import os, rasterio
	import numpy as np
	import multiprocessing as mp

	# list the data
	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	files = [ os.path.join( r, fn ) for r,s,files in os.walk( os.path.join( base_dir, 'CRU-TS40' ) ) 
					for fn in files if fn.endswith('.tif') and 'cru_ts40' in fn ]

	# rename
	pool = mp.Pool( 64 )
	out = pool.map( rename_cru_ts40, files )
	pool.close()
	pool.join()
