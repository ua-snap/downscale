# downscale the prepped cmip5 data used in running the TEM model (IEM)
# author: Michael Lindgren

if __name__ == '__main__':
	import glob, os, rasterio, itertools
	from functools import partial
	import downscale
	from downscale import preprocess
	import numpy as np
	import argparse

	# # parse the commandline arguments
	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-b", "--base_dir", action='store', dest='base_dir', type=str, help="base directory where data is stored in structured folders" )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="cmip5 model name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="cmip5 scenario name (exact)" )
	parser.add_argument( "-u", "--units", action='store', dest='units', type=str, help="cmip5 units name (exact)" )
	parser.add_argument( "-met", "--metric", action='store', dest='metric', type=str, help="cmip5 metric name (exact)" )
	parser.add_argument( "-lev", "--level", action='store', dest='level', type=int, help="optional level to extract for downscaling" )
	parser.add_argument( "-levn", "--level_name", action='store', dest='level_name', type=str, help="name of level variable" )
	
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

	if level is not None:
		level = float( level )

	# hardwired ARGS -- CMIP5
	project = 'ar5'
	interp = False
	find_bounds = False
	fix_clim = False
	aoi_mask = None # for precip data only
	anom = True # write out anoms (True) or not (False)
	
	# # # FOR TESTING # # # 
	# base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
	# variable = 'clt'
	# scenario = 'historical'
	# model = 'GFDL-CM3'
	# units = 'pct'
	# metric = 'mean'
	# level_name = None
	# level = None
	# # level = 1000 # mb / Pa
	# # level_name = 'plev'

	# # if level is not None:
	# # 	level = float( level )
	# # # # # # END TESTING # # # 
	
	# some setup args
	base_path = os.path.join( base_dir,'cmip5','prepped' )
	output_dir = os.path.join( base_dir, 'insolation_L48', 'downscaled_L48' )
	variables = [ variable ]
	scenarios = [ scenario ]
	models = [ model ]
	
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
		cru_cl20_varnames = {'hur':'reh', 'clt':'clt'} # we only support these variables for now...
		clim_path = os.path.join( base_dir, 'insolation_L48', 'climatologies', cru_cl20_varnames[variable] )
		filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
		filelist = [ i for i in filelist if '_14_' not in i ] # remove the GD ANNUAL _14_ file.
		baseline = downscale.Baseline( filelist )
		
		input_path = os.path.join( base_path, model, scenario, variable )
		output_path = os.path.join( output_dir, model, scenario, variable )

		if not os.path.exists( output_path ):
			os.makedirs( output_path )

		print( input_path )

		# list files for this set of downscaling -- one per folder
		fn, = glob.glob( os.path.join( input_path, '*.nc' ) )

		if 'historical' in scenario:
			historical = downscale.Dataset( fn, variable, model, scenario, project=project, units=units, 
							metric=metric, begin=1900, end=2005, level_name=level_name, level=level )
			future = None
		else:
			# get the historical data for anomalies
			historical_fn, = glob.glob( os.path.join( os.path.dirname( fn ).replace( scenario, 'historical' ), '*.nc' ) )
			historical = downscale.Dataset( historical_fn, variable, model, scenario, project=project, units=units, 
											metric=metric, begin=1900, end=2005, level_name=level_name, level=level )
			future = downscale.Dataset( fn, variable, model, scenario, project=project, units=units, metric=metric, 
											begin=2006, end=2100, level_name=level_name, level=level )

		# convert from Kelvin to Celcius
		if variable == 'tas':
			if historical:
				historical.ds[ variable ] = historical.ds[ variable ] - 273.15
				historical.ds[ variable ][ 'units' ] = units
			
			if future:
				future.ds[ variable ] = future.ds[ variable ] - 273.15
				future.ds[ variable ][ 'units' ] = units

		# DOWNSCALE
		mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
		clim_begin = '1961'
		clim_end = '1990'

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

		def round_data_clamp_hur( x ):
			x = round_data( x )
			x[ x < 0.0 ] = 0.0
			x[ x > 100.0 ] = 95.0 # per Stephanie McAfee
			return x

		def round_data_clamp_clt( x ):
			x = round_data( x )
			x[ x < 0.0 ] = 0.0
			x[ x > 100.0 ] = 100.0 # per Stephanie McAfee
			return x

		if variable == 'hur':
			post_downscale_function = round_data_clamp_hur
		elif variable == 'clt':
			post_downscale_function = round_data_clamp_clt
		else:
			post_downscale_function = round_data

		ar5 = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, future,
				downscaling_operation=downscaling_operation, mask=mask, mask_value=0, ncpus=64,
				src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
				post_downscale_function=post_downscale_function, varname=variable, modelname=modelname, 
				anom=anom, interp=interp, find_bounds=find_bounds, fix_clim=fix_clim, aoi_mask=aoi_mask )

		ar5.downscale( output_dir=output_path )
		