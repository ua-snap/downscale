# downscale the prepped cmip5 data downloaded using SYNDA for EPSCoR SC project
# author: Michael Lindgren -- June 09, 2016
def sort_files( files, split_on='_', elem_month=-2, elem_year=-1 ):
	'''
	sort a list of files properly using the month and year parsed
	from the filename.  This is useful with SNAP data since the standard
	is to name files like '<prefix>_MM_YYYY.tif'.  If sorted using base
	Pythons sort/sorted functions, things will be sorted by the first char
	of the month, which makes thing go 1, 11, ... which sucks for timeseries
	this sorts it properly following SNAP standards as the default settings.
	ARGUMENTS:
	----------
	files = [list] list of `str` pathnames to be sorted by month and year. usually from glob.glob.
	split_on = [str] `str` character to split the filename on.  default:'_', SNAP standard.
	elem_month = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
		default:-2. For SNAP standard.
	elem_year = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
		default:-1. For SNAP standard.
	RETURNS:
	--------
	sorted `list` by month and year ascending. 
	'''
	import pandas as pd
	months = [ int(fn.split('.')[0].split( split_on )[elem_month]) for fn in files ]
	years = [ int(fn.split('.')[0].split( split_on )[elem_year]) for fn in files ]
	df = pd.DataFrame( {'fn':files, 'month':months, 'year':years} )
	df_sorted = df.sort_values( ['year', 'month' ] )
	return df_sorted.fn.tolist()

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

if __name__ == '__main__':
	import glob, os, rasterio, itertools
	from functools import partial
	import downscale
	from downscale import preprocess
	import numpy as np
	import argparse

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-b", "--base_dir", action='store', dest='base_dir', type=str, help="base directory where data is stored in structured folders" )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="cmip5 model name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-mv", "--mean_variable", action='store', dest='mean_variable', type=str, help="cmip5 mean variable name (exact)" )
	parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="cmip5 scenario name (exact)" )
	parser.add_argument( "-u", "--units", action='store', dest='units', type=str, help="cmip5 units name (exact)" )
	parser.add_argument( "-met", "--metric", action='store', dest='metric', type=str, help="cmip5 metric name (exact)" )
	parser.add_argument( "-lev", "--level", action='store', dest='level', const=None, type=int, help="optional level to extract for downscaling" )
	parser.add_argument( "-levn", "--level_name", action='store', dest='level_name', const=None, type=str, help="name of level variable" )
	
	args = parser.parse_args()

	# unpack the args
	variable = args.variable
	scenario = args.scenario
	model = args.model
	units = args.units
	metric = args.metric
	base_dir = args.base_dir
	level = args.level
	level_name = args.level_name

	project = 'ar5'
	
	# # # # FOR TESTING # # # 
	# # base_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016'
	# base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data'
	# # fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/IPSL-CM5A-LR/rcp85/tasmax/tasmax_IPSL-CM5A-LR_rcp85_r1i1p1_2006_2100.nc'
	# # mean_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/IPSL-CM5A-LR/rcp85/tas/tas_IPSL-CM5A-LR_rcp85_r1i1p1_2006_2100.nc'
	# variable = 'tasmax'
	# mean_variable = 'tas'
	# scenario = 'rcp85'
	# model = 'IPSL-CM5A-LR'
	# units = 'C'
	# metric = 'mean'

	# some setup args
	base_path = os.path.join( base_dir,'cmip5','prepped' )
	# base_path = os.path.join( base_dir, 'downscaled' )
	downscaled_path = os.path.join( base_dir, 'downscaled' )
	raw_path = os.path.join( base_dir, 'cmip5', 'prepped' )
	output_dir = os.path.join( base_dir, 'downscaled_minmax_test' )
	variables = [ variable ]
	scenarios = [ scenario ]
	models = [ model ]
	anom = False # write out anoms (True) or not (False)

	# modelnames is simply the string name to put in the output filenaming if that differs from the modelname
	# used in querying the file which is the models list variable
	all_models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ] # temp for distributed run
	modelnames = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4' ]

	modelnames = dict( zip( all_models, modelnames ) )
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )

	os.chdir( output_dir )

	for variable, model, scenario in itertools.product( variables, models, scenarios ):
		modelname = modelnames[ model ]
		# SETUP BASELINE
		clim_path = os.path.join( base_dir, 'downscaled', model, scenario, mean_variable )
		filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
		# sort these files
		filelist = only_years( sort_files( filelist ), begin=2006, end=2100 )

		# filelist = [ i for i in filelist if '_14_' not in i ] # remove the GD ANNUAL _14_ file.
		baseline = downscale.Baseline( filelist )
		
		input_path = os.path.join( base_path, model, scenario, variable )
		output_path = os.path.join( output_dir, model, scenario, variable )

		if not os.path.exists( output_path ):
			os.makedirs( output_path )

		print( input_path )

		# # NOTE:
		# ALL DATA FALL INTO THE historical ARGUMENT WITH THIS DOWNSCALING DUE TO NOT USING A CLIMATOLOGY TO 
		# GENERATE ANOMALIES WE ARE GENERATING DELTAS OF MIN/MAX FROM MEAN OF THE SAME VARIABLE GROUP ie. tas/tasmax/tasmin
		# list files for this set of downscaling -- one per folder
		if scenario == 'historical':
			begin = 1900
			end = 2005
		else:
			begin = 2006
			end = 2100

		fn, = glob.glob( os.path.join( input_path, '*.nc' ) )
		historical = downscale.Dataset( fn, variable, model, scenario, project=project, units=units, metric=metric, begin=begin, end=end )
		# mean data -- hacky...
		mean_fn, = glob.glob( os.path.join( input_path.replace( variable, mean_variable ), '*.nc' ) )
		mean_ds = downscale.Dataset( mean_fn, mean_variable, model, scenario, project=project, units=units, metric=metric, begin=begin, end=end )
		future = None

		# DOWNSCALE
		mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
		
		# these have absolutely no effect but are here since they are a required variable to the super class DeltaDownscale...
		# we need a way to make this more nimble as this is not ideal...
		clim_begin = '1961'
		clim_end = '1990'

		# rounding switch
		if variable == 'pr':
			rounder = np.rint
			downscaling_operation = 'mult'
		elif variable in ['hur','cld','clt']:
			rounder = partial( np.round, decimals=1 )
			downscaling_operation = 'mult'
		else:
			rounder = partial( np.round, decimals=1 )
			downscaling_operation = 'add'

		def round_it( x, mask ):
			arr = np.ma.masked_array( data=x, mask=mask )
			return rounder( arr )

		round_data = partial( round_it, mask=( mask==0 ) )

		def round_data_clamp( x ):
			x[ x < 0.0 ] = 0.0
			x[ x > 100.0 ] = 100.0
			return round_data( x )

		if variable == 'hur' or variable == 'clt':
			post_downscale_function = round_data_clamp
		else:
			post_downscale_function = round_data

		ar5 = downscale.DeltaDownscaleMinMax( baseline=baseline, clim_begin=clim_begin, clim_end=clim_end, historical=historical, future=future,
					downscaling_operation=downscaling_operation, mask=mask, mask_value=0, ncpus=32,
					src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
					post_downscale_function=post_downscale_function, varname=variable, 
					modelname=modelname, anom=anom, mean_ds=mean_ds, mean_variable=mean_variable )

		ar5.downscale( output_dir=output_path )
		