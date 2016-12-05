# data splitter dof/dot/logs
def read_split( fn ):
	'''
	read and divide up the 3 band dof/dot/logs data
	into individual variable rasters
	'''
	import os
	import numpy as np

	dirname, basename = os.path.split( fn )

	with rasterio.open( fn ) as rst:
		meta = rst.meta
		meta.update( compress='lzw', count=1 )
		meta.pop( 'transform' )

		varnames = ['dof', 'dot', 'logs'] # [hardwired] order is important here.
		for variable, arr in zip( varnames, rst.read() ):
			out_dirname = os.path.join( dirname.replace( '/dot_dof_grow', '' ), variable ) # [hardwired]
			out_basename = basename.replace( 'dot_dof_grow', variable )
			output_filename = os.path.join( out_dirname, out_basename )

			try:
				if not os.path.exists( out_dirname ):
					os.makedirs( out_dirname )
			except:
				pass

			with rasterio.open( output_filename, 'w', **meta ) as out:
				out.write( arr, 1 )
	return 1

if __name__ == '__main__':
	import os, rasterio
	import numpy as np
	from pathos.mp_map import mp_map

	# base path 
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/decadal_monthlies'

	# list the data we want to split up
	files = [ os.path.join( r, f ) for r,s,files in os.walk( base_path ) for f in files if f.endswith( '.tif' ) and 'dot_dof_grow' in f ]

	# split em up
	done = mp_map( read_split, files, nproc=32 )
