# downscale the prepped cmip5 data downloaded using SYNDA for EPSCoR SE project
# author: Michael Lindgren -- June 09, 2016
if __name__ == '__main__':
	import glob, os, rasterio, itertools
	import downscale
	from downscale import preprocess
	import argparse

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="cmip5 model name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="cmip5 scenario name (exact)" )	
	args = parser.parse_args()

	# unpack the args
	variable = args.variable
	scenario = args.scenario
	model = args.model

	project = 'ar5'
	units = 'C'
	metric = 'mean'
	
	# # # # # TESTING # # # 
	# variable = 'tasmin'
	# scenario = 'rcp26'
	# model = 'CCSM4'

	# some setup args
	base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prepped_cmip5'
	output_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_v2'
	variables = [ variable ] # [ 'tasmin', 'tasmax' ]
	scenarios = [ scenario ] # [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	models = [ model ] # , 'MRI-CGCM3''GISS-E2-R', 'GFDL-CM3', 'CCSM4'
	# [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ]

	# modelnames is simply the string name to put in the output filenaming if that differs from the modelname
	# used in querying the file which is the models list variable
	all_models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ] # temp for distributed run
	modelnames = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4' ]

	modelnames = dict(zip( all_models, modelnames ))
 	# # #

	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )

	# open a log file to find out where we are messing up
	log = open( os.path.join( output_dir, 'log_file_downscale.txt' ), 'w' )

	for variable, model, scenario in itertools.product( variables, models, scenarios ):
		modelname = modelnames[ model ]
		# SETUP BASELINE
		if variable == 'tasmin':
			v = 'tmin'
		elif variable == 'tasmax':
			v = 'tmax'
		else:
			NotImplementedError( 'only tasmin and tasmax are currently supported' )

		clim_path = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/akcan', v )
		filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
		baseline = downscale.Baseline( filelist )

		input_path = os.path.join( base_dir, model, scenario, variable )
		output_path = os.path.join( output_dir, model, scenario, variable )

		if not os.path.exists( output_path ):
			os.makedirs( output_path )

		print( input_path )

		# list files for this set of downscaling -- one per folder
		fn, = glob.glob( os.path.join( input_path, '*.nc' ) )

		if 'historical' in scenario:
			historical = downscale.Dataset( fn, variable, model, scenario, project=project, units=units, metric=metric )
			future = None # no need for futures here....
		else:
			# get the historical data for anomalies
			historical_fn, = glob.glob( os.path.join( os.path.dirname( fn ).replace( scenario, 'historical' ), '*.nc' ) )
			historical = downscale.Dataset( historical_fn, variable, model, scenario, project=project, units=units, metric=metric )
			future = downscale.Dataset( fn, variable, model, scenario, project=project, units=units, metric=metric )

		# convert from Kelvin to Celcius
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

		ar5 = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, future, \
				downscaling_operation='add', mask=mask, mask_value=0, ncpus=32, \
				src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
				post_downscale_function=None, varname=variable, modelname=modelname ) # -9999.0

		ar5.downscale( output_dir=output_path )
		
