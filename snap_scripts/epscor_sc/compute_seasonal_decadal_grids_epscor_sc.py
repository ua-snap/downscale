# def sort_files( files, split_on='_', elem_month=-2, elem_year=-1 ):
# 	'''
# 	sort a list of files properly using the month and year parsed
# 	from the filename.  This is useful with SNAP data since the standard
# 	is to name files like '<prefix>_MM_YYYY.tif'.  If sorted using base
# 	Pythons sort/sorted functions, things will be sorted by the first char
# 	of the month, which makes thing go 1, 11, ... which sucks for timeseries
# 	this sorts it properly following SNAP standards as the default settings.

# 	ARGUMENTS:
# 	----------
# 	files = [list] list of `str` pathnames to be sorted by month and year. usually from glob.glob.
# 	split_on = [str] `str` character to split the filename on.  default:'_', SNAP standard.
# 	elem_month = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
# 		default:-2. For SNAP standard.
# 	elem_year = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
# 		default:-1. For SNAP standard.

# 	RETURNS:
# 	--------
# 	sorted `list` by month and year ascending. 

# 	'''
# 	import pandas as pd
# 	months = [ int(fn.split('.')[0].split( split_on )[elem_month]) for fn in files ]
# 	years = [ int(fn.split('.')[0].split( split_on )[elem_year]) for fn in files ]
# 	df = pd.DataFrame( {'fn':files, 'month':months, 'year':years} )
# 	df_sorted = df.sort_values( ['year', 'month' ] )
# 	return df_sorted.fn.tolist()

def only_years( files, begin=1901, end=2100, split_on='_', elem_year=-1 ):
	'''
	return new list of filenames where they are truncated to begin:end

	ARGUMENTS:
	----------
	files = [list] list of `str` pathnames to be sorted by month and year. usually from glob.glob.
	begin = [int] four digit integer year of the begin time default:1901
	end = [int] four digit integer year of the end time default:2100
	split_on = [str] `str` character to split the filename on.  default:'_', SNAP standard.
	elem_year = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
		default:-1. For SNAP standard.

	RETURNS:
	--------
	sliced `list` to begin and end year.
	'''
	import pandas as pd
	years = [ int(fn.split('.')[0].split( split_on )[elem_year]) for fn in files ]
	df = pd.DataFrame( { 'fn':files, 'year':years } )
	df_slice = df[ (df.year >= begin ) & (df.year <= end ) ]
	return df_slice.fn.tolist()

# def get_month_seaon( fn ):
# 	# seasons
# 	seasonal_lookup = { 1:'DJF', 2:'DJF', 3:'MAM', 4:'MAM', 5:'MAM', \
# 						6:'JJA', 7:'JJA', 8:'JJA',\
# 						 9:'SON', 10:'SON', 11:'SON', 12:'DJF' }

# 	fn = os.path.basename( fn )
# 	month, year = fn.replace( '.tif', '' ).split( '_' )[-2:]
# 	return seasonal_lookup[ int(month) ]

def get_year( fn ):
	fn = os.path.basename( fn )
	month, year = fn.replace( '.tif', '' ).split( '_' )[-2:]
	return year

def get_season( fn ):
	fn = os.path.basename( fn )
	season, year = fn.replace( '.tif', '' ).split( '_' )[-2:]
	return season

def read_raster( fn, band=1 ):
	'''
	clean way to open / read and properly close a GTiff
	'''
	import rasterio
	with rasterio.open( fn ) as out:
		arr = out.read( band )
	return arr

def calc_decadal_from_annual( season_name, files, output_path, agg_metric, *args, **kwargs ):
	'''
	calculate seasonal means
	'''
	from functools import partial

	years = [ int( get_year( fn ) ) for fn in files ]
	year = str( max( years ) )
	fn = files[0]
	with rasterio.open( fn ) as rst:
		mask = rst.read_masks( 1 )
		meta = rst.meta
		
	if 'transform' in meta.keys():
		meta.pop( 'transform' )

	meta.update( compress='lzw' )

	# NOTE HOW FOR DECADALS WE ARE ONLY PERFORMING MEANS!!!!  THIS IS IMPORTANT!!!
	metric_switch = { 'mean':np.mean, 'total':np.mean, 'min':np.min, 'max':np.max }

	variable, metric, units, project, model, scenario = os.path.basename( fn ).split( '.' )[0].split( '_' )[:-2]
	arr = np.array([ read_raster( i ) for i in files ])
	arr = metric_switch[ agg_metric ]( arr, axis=0 )

	# round the data 
	if variable == 'pr':
		# truncate to whole number
		rounder = partial( np.round, decimals=0 )
	else:
		# round to 2 decimals
		rounder = partial( np.round, decimals=1 )

	def round_it( x, mask ):
		arr = np.ma.masked_array( x, mask == 0 )
		return rounder( arr )

	# round the data -- though for > 0 it does nothing
	round_data = partial( round_it, mask=mask )
	arr = round_data( arr ).data
	arr[ mask == 0 ] = meta[ 'nodata' ]

	decade_out = str(year)[:3] + '0s'
	output_filename = os.path.join( output_path, model, scenario, variable, '_'.join([variable, agg_metric, units, project, model, scenario, season_name, decade_out]) + '.tif' )

	dirname = os.path.dirname( output_filename )
	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( arr, 1 )

	return output_filename

def make_decadal_seasonal( base_path, output_path, variable, model, scenario, decade, ncpus, agg_metric ):
	'''
	function to calculate and output mean seasonal monthly data across decades
	
	ARGUMENTS:
	----------
	base_path = [  ]  
	output_path = [  ]  
	model = [  ]  
	scenario = [  ]  
	variable = [  ]  
	begin = [  ]  
	end = [  ]  
	ncpus = [  ]  

	RETURNS
	-------
	output_directory of newly produced GeoTiffs if successful. else, error.

	'''
	decade_begin, decade_end = decade

	# modeled data
	files = glob.glob( os.path.join( base_path, model, scenario, variable, '*' + agg_metric + '*.tif' ) )
	files = only_years( files, begin=decade_begin, end=decade_end, split_on='_', elem_year=-1 )

	# season_names = [ get_month_seaon( fn ) for fn in files ]
	years = [ int(get_year( fn )) for fn in files ]

	# min / max years
	start_year =  str( min(years) )
	end_year = str( max(years) )

	seasons = [ get_season( fn ) for fn in files ]

	# drop data for start_year JF and end_year this is useful for annuals, but not really decadals
	# files = [ fn for fn in files if not '_'.join([ '01',start_year ]) in fn if not '_'.join([ '02',start_year ]) in fn if not '_'.join([ '12',end_year ]) in fn ]
	files = pd.Series( files )

	grouped_seasons = files.groupby( seasons )

	args = [ ( season_name, file_group.tolist(), output_path, agg_metric ) for season_name, file_group in grouped_seasons ]

	_ = mp_map( wrap, args, nproc=ncpus )
	return args

def wrap( x ):
	''' 
	multiprocessing wrapper for clean 
	argument handling without lambda 
	'''
	return calc_decadal_from_annual( *x )

if __name__ == '__main__':
	import os, glob, itertools, rasterio
	import xarray as xr
	import pandas as pd
	import numpy as np
	from pathos.mp_map import mp_map
	import argparse

	'''
	this tool assumes that the data are stored in a directory structure as follows:
	
	base_path
		model
			scenario
				variable
					FILES
	'''

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-b", "--base_path", action='store', dest='base_path', type=str, help="path to the directory where the downscaled modeled data are stored" )
	parser.add_argument( "-o", "--output_path", action='store', dest='output_path', type=str, help="path to the output directory" )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="model name (exact)" )
	parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="scenario name (exact)" )
	parser.add_argument( "-p", "--project", action='store', dest='project', type=str, help="project name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-am", "--agg_metric", action='store', dest='agg_metric', type=str, help="string name of the metric to compute the decadal summary - mean, max, min, total" )
	parser.add_argument( "-nc", "--ncpus", action='store', dest='ncpus', type=int, help="number of cpus to use in multiprocessing" )	
	args = parser.parse_args()

	# unpack for cleaner var access:
	base_path = args.base_path
	output_path = args.output_path
	model = args.model
	scenario = args.scenario
	project = args.project
	variable = args.variable
	ncpus = args.ncpus
	agg_metric = args.agg_metric

	# # # # # FOR TESTING
	# base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_FINAL_OCT_TESTING'
	# output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_FINAL_OCT_TESTING'
	# model = 'GFDL-CM3'
	# scenario = 'rcp60'
	# project = 'cmip5'
	# variable = 'pr'
	# agg_metric = 'total'
	# ncpus = 32

	# # # # # # # # # # #

	# switches to deal with different date groups.  Hardwired to CMIP5 and CRU TS323 currently.
	cmip_switch = { 'historical':(1900,2005), 'rcp26':(2006,2100), 'rcp45':(2006,2100), 'rcp60':(2006,2100), 'rcp85':(2006,2100) }
	cru_switch = { 'historical':(1901,2014) }
	project_switch = { 'cmip5':cmip_switch, 'cru':cru_switch }

	begin, end = project_switch[ project ][ scenario ]
	decades = list( sorted( set( [ (int(str(i)[:3]+'0'), int(str(i)[:3]+'9')) for i in range( begin, end ) ] ) ) )

	for decade in decades:
		print( 'running: {} {} {} {}'.format( model, variable, scenario, decade ) )
		done = make_decadal_seasonal( base_path, output_path, variable, model, scenario, decade, ncpus, agg_metric )





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # FOR TESTING THEN REMOVE
# # # setup args
# # cmip5
# import subprocess, os
# # base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_clipped'
# base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/annual_seasonals'
# output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/decadal_seasonals'
# ncpus = 32
# # project = 'cru'
# project = 'cmip5'
# variables = [ 'tasmin', 'tasmax', 'tas', 'pr' ]
# # models = [ 'ts323' ]
# models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4', '5ModelAvg' ]
# # scenarios = [ 'historical']
# scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]

# for model in models:
# 	for scenario in scenarios:
# 		for variable in variables:
# 			agg_metric = 'mean'

# 			# if variable == 'pr':
# 			# 	agg_metric = 'total'
# 			# else:
# 			# 	agg_metric = 'mean'
# 			os.chdir( '/workspace/UA/malindgren/repos/downscale/snap_scripts' )
# 			command = ' '.join([ 'ipython', 'compute_seasonal_decadal_grids_epscor_se.py', '--', '-b', base_path, '-o ', output_path, '-m ', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])

# 			os.system( command )

# import subprocess, os
# # base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_clipped'
# base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_clipped'
# output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_outputs/decadal_seasonal'
# ncpus = 32
# # project = 'cru'
# project = 'cmip5'
# variables = [ 'tasmin', 'tasmax', 'tas', 'pr' ]
# # models = [ 'ts323' ]
# models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4', '5ModelAvg' ]
# # scenarios = [ 'historical']
# scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]

# for model in models:
# 	for scenario in scenarios:
# 		for variable in variables:
# 			# agg_metric = 'mean'
# 			if variable == 'pr':
# 				agg_metric = 'total'
# 			else:
# 				agg_metric = 'mean'
# 			os.chdir( '/workspace/UA/malindgren/repos/downscale/snap_scripts' )
# 			command = ' '.join([ 'ipython', 'compute_seasonal_decadal_grids_epscor_se.py', '--', '-b', base_path, '-o ', output_path, '-m ', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
# 			os.system( command )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # cru
# import subprocess, os
# base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_minmax'
# # base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_clipped'
# output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_outputs_minmax/decadal_seasonal'
# ncpus = 32
# project = 'cru'
# # project = 'cmip5'
# variables = ['pr'] #[ 'tasmin', 'tasmax', 'tas', 'pr' ]
# models = [ 'ts323' ]
# # models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4', '5ModelAvg' ]
# scenarios = [ 'historical']
# # scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]

# for model in models:
# 	for scenario in scenarios:
# 		for variable in variables:
# 			agg_metric = 'mean'
# 			# if variable == 'pr':
# 			# 	agg_metric = 'total'
# 			# else:
# 			# 	agg_metric = 'mean'
# 			os.chdir( '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc' )
# 			command = ' '.join([ 'ipython', 'compute_seasonal_decadal_grids_epscor_sc.py', '--', '-b', base_path, '-o ', output_path, '-m ', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
# 			os.system( command )

# import subprocess, os
# base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_clipped'
# # base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_clipped'
# output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_outputs/decadal_seasonal'
# ncpus = 32
# project = 'cru'
# # project = 'cmip5'
# variables = [ 'tasmin', 'tasmax', 'tas', 'pr' ]
# models = [ 'ts323' ]
# # models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4', '5ModelAvg' ]
# scenarios = [ 'historical']
# # scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]

# for model in models:
# 	for scenario in scenarios:
# 		for variable in variables:
# 			# agg_metric = 'mean'
# 			if variable == 'pr':
# 				agg_metric = 'total'
# 			else:
# 				agg_metric = 'mean'
# 			os.chdir( '/workspace/UA/malindgren/repos/downscale/snap_scripts' )
# 			command = ' '.join([ 'ipython', 'compute_seasonal_decadal_grids_epscor_se.py', '--', '-b', base_path, '-o ', output_path, '-m ', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
# 			os.system( command )





# # # # # OLDER REMOVE
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # FOR TESTING THEN REMOVE
# # # setup args
# # import subprocess, os
# # base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_clipped'
# # # base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_clipped'
# # output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_outputs/monthly_decadal_seasonal'
# # ncpus = 32
# # project = 'cru'
# # # project = 'cmip5'
# # variables = [ 'tasmin', 'tasmax', 'tas', 'pr' ]
# # models = [ 'ts323' ]
# # # models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4', '5ModelAvg' ]
# # scenarios = [ 'historical']
# # # scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]

# # for model in models:
# # 	for scenario in scenarios:
# # 		for variable in variables:
# # 			agg_metric = 'mean'
# # 			if variable == 'pr':
# # 				agg_metric = 'total'
# # 			else:
# # 				agg_metric = 'mean'
# # 			os.chdir( '/workspace/UA/malindgren/repos/downscale/snap_scripts' )
# # 			command = ' '.join([ 'python','compute_seasonal_decadal_grids_epscor_se.py', '-b', base_path, '-o ', output_path, '-m ', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
# # 			os.system( command )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

