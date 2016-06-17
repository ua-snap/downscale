# Generate 5 Model Averages for EPSCoR SE project
# author: Michael Lindgren -- June 09, 2016

def list_files( input_path, begin=1900, end=2005 ):
	''' list the files in the range of begin, end '''
	files = sorted( glob.glob( os.path.join( input_path, '*.tif' ) ) )
	# var, units, metric, model, month, year = os.path.basename( fn ).split( '.' )[0].split( '_' )
	return [ fn for fn in files if int(os.path.basename( fn ).split( '.' )[0].split( '_' )[-1]) in range( begin, end+1 ) ]

def generate( files ):
	rst = rasterio.open( files[0] )
	meta = rst.meta
	meta.update( compress='lzw' )
	mask = rst.read_masks( 1 )

	if 'transform' in meta.keys():
		meta.pop( 'transform' )

	arr = np.mean([ rasterio.open( fn ).read( 1 ) for fn in files ], axis=0 )
	arr[ mask == 0 ] = meta[ 'nodata' ]

	fn, = [ i for i in files if 'CCSM4' in i ]
	output_filename = fn.replace( 'CCSM4', '5ModelAvg' )
	new_path = os.path.dirname( output_filename )
	# new directory -- if needed:
	if not os.path.exists( new_path ):
		os.makedirs( new_path )

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( arr, 1 )
	return output_filename

if __name__ == '__main__':
	import glob, os, rasterio, itertools
	import numpy as np
	from pathos import multiprocessing

	# some setup args
	base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5'
	variables = [ 'tasmin', 'tasmax' ]
	scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ]

	for variable, scenario in itertools.product( variables, scenarios ):
		if scenario == 'historical':
			begin = 1900
			end = 2005
		else:
			begin = 2006
			end = 2100

		# list the files we want
		input_files = [ list_files( os.path.join( base_dir, model, scenario, variable ), begin, end ) for model in models ]
		grouped = zip( *input_files )

		pool = multiprocessing.Pool( 32 )
		output_filenames = pool.map( generate, grouped )
		pool.close()
		pool.join()
