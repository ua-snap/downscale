# # # # # # # # # # # # # # # # # # # # # # # # 
# Generate 5 Model Averages --> 10 min data
# author: Michael Lindgren -- SEPT 2017
# # # # # # # # # # # # # # # # # # # # # # # # 

def list_files( input_path, begin=1900, end=2005 ):
	''' list the files in the range of begin, end '''
	files = sorted( glob.glob( os.path.join( input_path, '*.tif' ) ) )
	# var, units, metric, model, month, year = os.path.basename( fn ).split( '.' )[0].split( '_' )
	return [ fn for fn in files if int(os.path.basename( fn ).split( '.' )[0].split( '_' )[-1]) in range( begin, end+1 ) ]

def generate( files ):
	from functools import partial
	# from pathos.mp_map import mp_map

	rst = rasterio.open( files[0] )
	meta = rst.meta
	meta.update( compress='lzw' )
	mask = rst.read_masks( 1 )

	if 'transform' in meta.keys():
		meta.pop( 'transform' )

	# arr_group = mp_map( lambda x:rasterio.open( x ).read( 1 ), files, nproc=32 )
	arr_group = np.array([ rasterio.open( fn ).read( 1 ) for fn in files ])
	group_mean = np.mean( arr_group, axis=0 )
	group_mean[ mask == 0 ] = rst.nodata # add back the oob -3.39999995e+38 

	fn, = [ i for i in files if 'NCAR-CCSM4' in i ]
	output_filename = fn.replace( 'NCAR-CCSM4', '5ModelAvg' )
	new_path = os.path.dirname( output_filename )
	try:
		# new directory -- if needed:
		if not os.path.exists( new_path ):
			os.makedirs( new_path )
	except:
		pass

	# round the data 
	# dynamically get the first element of the filename which is variable
	variable = os.path.basename(files[0]).split( '_' )[0]
	
	if variable == 'pr':
		# truncate to whole number
		rounder = np.rint
	else:
		# round to 2 decimals
		rounder = partial( np.round, decimals=1 )

	def round_it( x, mask ):
		arr = np.ma.masked_array( x, mask == 0 )
		return rounder( arr )

	# round the data -- though for > 0 it does nothing
	group_mean = round_it( group_mean, mask=mask )
	group_mean.fill_value = rst.nodata

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( group_mean.filled(), 1 )
	return output_filename

if __name__ == '__main__':
	import glob, os, rasterio, itertools
	import numpy as np
	from pathos.mp_map import mp_map
	import argparse

	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-b", "--base_dir", action='store', dest='base_dir', type=str, help="base directory where data is stored in structured folders" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="scenario name (exact)" )
	args = parser.parse_args()
	
	# unpack
	variable = args.variable
	base_dir = args.base_dir
	scenario = args.scenario

	# some setup args
	base_dir = os.path.join( base_dir, 'downscaled_10min' )
	variables = [ variable ]
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4' ]

	for variable in variables:
		if scenario == 'historical':
			begin = 1860
			end = 2005
		else:
			begin = 2006
			end = 2100

		# list the files we want
		input_files = [ list_files( os.path.join( base_dir, model, scenario, variable ), begin, end ) for model in models ]
		grouped = zip( *input_files )

		# run it in parallel
		output_filenames = mp_map( generate, grouped, nproc=32 )
