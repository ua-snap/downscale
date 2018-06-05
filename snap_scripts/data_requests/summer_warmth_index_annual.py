# # # # # # # 
# CALCULATE ANNUAL Summer Warmth Index from SNAP standard 2km AKCAN outputs
# # # # # # # 

def list_files( base_dir, wildcard='tas_*.tif' ):
	import glob
	out = []
	for root, subs, files in os.walk( base_dir ):
		if len( [i for i in files if 'tas_' in i ] ) > 0:
			out = out + glob.glob(os.path.join( root, wildcard ) )
	return out

def make_grouper( x ):
	group = os.path.basename( x ).split('.')[0].split('_')
	year = group[-1]
	return year

def swi( files, output_path ):
	import rasterio
	import numpy as np

	months = [ fn.split('_')[-2] for fn in files ]
	files = [ fn for fn, month in zip( files, months ) if month in ['05', '06', '07', '08', '09'] ]

	arr = np.array( [ rasterio.open(i).read(1) for i in files ] )
	arr[ arr < 0 ] = 0
	arr = np.sum( arr, axis=0 )
	
	variable, metric, units, project, model, \
		scenario, month, year = \
			os.path.basename( files[0] ).split( '.' )[0].split( '_' )

	# make sure the year is a decadal year for grouping if annual...
	# this works with both decadal and annual monthlies.
	rst = rasterio.open( files[0] )
	mask = rst.read_masks( 1 )
	meta = rst.meta

	meta.update( compress='lzw' )
	output_filename = os.path.join( output_path, model, scenario, 'swi', '_'.join(['swi', 'cumulative', units, project, model, scenario, year])+'.tif' )

	try:
		dirname = os.path.dirname( output_filename )
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass

	with rasterio.open( output_filename, 'w', **meta ) as out:
		arr[ mask == 0 ] = meta[ 'nodata' ]
		out.write( np.round( arr, 1 ), 1 )
	return output_filename

if __name__ == '__main__':
	import os
	import pandas as pd
	import multiprocessing as mp
	from functools import partial
	import argparse

	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-b", "--base_path", action='store', dest='base_path', type=str, help="path to the directory where the downscaled modeled data are stored" )
	parser.add_argument( "-o", "--output_path", action='store', dest='output_path', type=str, help="path to the output directory" )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="model name (exact)" )
	parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="scenario name (exact)" )
	parser.add_argument( "-p", "--project", action='store', dest='project', type=str, help="project name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )

	# unpack args for cleaner access
	args = parser.parse_args()
	base_path = args.base_path
	output_path = args.output_path
	model = args.model
	scenario = args.scenario
	project = args.project
	variable = 'tas' # swi only on tas data

	# list data
	l = pd.Series( list_files( os.path.join( base_path, model, scenario, variable ) ) )

	# group it by year
	groups = l.apply( make_grouper )
	grouped = l.groupby( groups )
	data_grouped = [ sorted(group) for name, group in grouped ]

	# parallel
	swi_f = partial( swi, output_path=output_path )
	pool = mp.Pool(32)
	done = pool.map( swi_f, data_grouped )
	pool.close()
	pool.join()


# # # # # NOTES:
	# # TESTING STUFF
	# base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	# model = 'GFDL-CM3'
	# scenario = 'rcp60'
	# variable = 'tas'
	# begin = 2006
	# end = 2100
	# ncpus = 32
	# metric = 'mean'
	# project = 'ar5'
	# output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/derived'
	# # # # # # # # # # 

# # run example
# import os, itertools, subprocess
# import 
# base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
# output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/derived'
# models = [ 'GFDL-CM3', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'NCAR-CCSM4', '5ModelAvg' ]
# scenarios = [ 'historical', 'rcp45', 'rcp60', 'rcp85' ]
# project = 'ar5'
# for model, scenario in itertools.product( models, scenarios):
# 	_ = subprocess.call(['summer_warmth_index_annual.py', '-b', base_path, '-o', output_path, '-m', model, '-s', scenario, '-p', project, '-v', 'tas' ])

