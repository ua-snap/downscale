def get_monyear( fn ):
	'''
	specialized function to split the filename 
	pattern of SNAP data and extract 
	year and month information from it.
	'''
	fn, ext = os.path.splitext( fn )
	return fn.split( '_' )[-2:]

def get_decade( fn ):
	''' get decade from snap-style filename '''
	month, year = get_monyear( fn )
	return year[:3]+'0s'

def get_month_seaon( fn ):
	''' method to use a month number as a way to determine seasons '''
	# seasons
	seasonal_lookup = { 1:'DJF', 2:'DJF', 3:'MAM', 4:'MAM', 5:'MAM', \
						6:'JJA', 7:'JJA', 8:'JJA',\
						 9:'SON', 10:'SON', 11:'SON', 12:'DJF' }

	fn = os.path.basename( fn )
	month, year = get_monyear( fn )
	return seasonal_lookup[ int(month) ]

def read_raster( fn, band=1 ):
	'''
	clean way to open / read and properly close a GTiff
	'''
	import rasterio
	with rasterio.open( fn ) as out:
		arr = out.read( band )
	return arr

def write_raster( arr, out_fn, meta ):
	import os
	meta.update( compress='lzw', dtype='float32' )
	dirname, basename = os.path.split( out_fn )
	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass

	with rasterio.open( out_fn, 'w', **meta ) as out:
		out.write( arr.astype( np.float32 ), 1)
	return out_fn

def aggregate( files_df, metric ):
	''' aggregate to a np method metric passed '''
	import numpy as np
	files = files_df.files.tolist()
	decade = files_df[ 'decades' ].iloc[0]
	lookup = { 'mean':np.mean, 'total':np.sum }
	arr = lookup[ metric ]([ rasterio.open(i).read(1) for i in files ], axis=0)
	return ( decade, files, arr )


if __name__ == '__main__':
	import glob, os, rasterio
	import pandas as pd
	import numpy as np
	from pathos.mp_map import mp_map
	from functools import partial
	import argparse

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-b", "--base_dir", action='store', dest='base_dir', type=str, help="path to the directory where the downscaled modeled data are stored" )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="model name (exact)" )
	parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="scenario name (exact)" )
	parser.add_argument( "-p", "--project", action='store', dest='project', const=None, type=str, help="project name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-am", "--agg_metric", action='store', dest='agg_metric', type=str, help="string name of the metric to compute the decadal summary - mean, max, min, total" )
	parser.add_argument( "-nc", "--ncpus", action='store', dest='ncpus', type=int, help="number of cpus to use in multiprocessing" )	
	args = parser.parse_args()

	# unpack for cleaner var access:
	base_dir = args.base_dir
	model = args.model
	scenario = args.scenario
	project = args.project
	variable = args.variable
	metric = args.agg_metric
	ncpus = args.ncpus
	
	files = glob.glob( os.path.join( base_dir, 'derived_grids','annual_seasonals', 
						model, scenario, variable, '*.tif' ) )

	# split to decades
	decades = [ get_decade( fn ) for fn in files ]
	seasons = [ get_monyear( fn )[0] for fn in files ] # season name in month spot
	df = pd.DataFrame( {'decades':decades, 'files':files, 'seasons':seasons})
	file_groups = [ sub_df for idx,sub_df in df.groupby( [decades, seasons] ) if len(sub_df.files) == 10 ]
	
	# compute
	_aggregate = partial( aggregate, metric=metric )
	agg_groups = mp_map( _aggregate, file_groups, nproc=ncpus )

	# make some output metadata
	meta = rasterio.open( files[0] ).meta
	meta.update( compress='lzw' )
	if 'transform' in meta.keys():
		meta.pop( 'transform' )

	mask = rasterio.open( files[0] ).read_masks( 1 )

	# make output filenames
	args = []
	for decade, filenames, arr in agg_groups:
		dirname, basename = os.path.split( filenames[0] )
		basename, ext = os.path.splitext( basename )
		variable, metric, units, project, model, scenario, season, year = basename.split( '_' )

		if variable == 'pr' and metric == 'mean':
			# handle pr mean_total
			metric = 'mean_total'

		new_basename = '_'.join([ variable, metric, units, project, model, scenario, decade ]) + ext

		if metric == 'mean_total':
			# handle pr mean_total
			metric = 'mean'		
		
		out_fn = os.path.join( dirname.replace( 'annual_seasonals', 'decadal_annual_seasonals' ), new_basename )
		arr = np.rint( arr )
		arr[ mask == 0 ] = meta['nodata']
		args = args + [( arr, out_fn, meta )]

	# write it out
	final = mp_map( lambda x: write_raster( *x ), args, nproc=32, **{'meta':meta} )


# # TESTING STUFF
# base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data'
# model = 'GFDL-CM3'
# scenario = 'rcp60'
# variable = 'pr'
# begin = 2006
# end = 2100
# ncpus = 32
# metric = 'mean'
# project = 'ar5'

