# downscale the prepped cmip5 data downloaded using SYNDA for EPSCoR SC project
# author: Michael Lindgren -- June 09, 2016 (UPDATED: September 21, 2016 -- [ML])
if __name__ == '__main__':
	import glob, os, rasterio, itertools
	from functools import partial
	import downscale
	from downscale import preprocess, Mask, utils
	import argparse
	import numpy as np

	# # parse the commandline arguments
	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-b", "--base_dir", action='store', dest='base_dir', type=str, help="base directory where data is stored in structured folders" )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="cmip5 model name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="cmip5 scenario name (exact)" )
	parser.add_argument( "-u", "--units", action='store', dest='units', type=str, help="cmip5 units name (exact)" )
	parser.add_argument( "-met", "--metric", action='store', dest='metric', type=str, help="cmip5 metric name (exact)" )
	args = parser.parse_args()

	# unpack the args
	variable = args.variable
	scenario = args.scenario
	model = args.model
	units = args.units
	metric = args.metric
	base_dir = args.base_dir

	# AOI MASK -- HARDWIRE -- PCLL for CMIP5
	aoi_mask_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/akcan_template/akcan_aoi_mask_PCLL.shp'
	project = 'ar5'
	
	# # # # FOR TESTING # # # 
	# base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data'
	# variable = 'pr'
	# scenario = 'rcp60'
	# model = 'GFDL-CM3'
	# units = 'mm'
	# metric = 'total'

	# some setup args
	base_path = os.path.join( base_dir,'cmip5','prepped' )
	output_dir = os.path.join( base_dir, 'downscaled' )
	variables = [ variable ]
	scenarios = [ scenario ]
	models = [ model ]
	anom = False # write out anoms (True) or not (False)
	interp = False # interpolate across space -- Low Res
	find_bounds = False

	# modelnames is simply the string name to put in the output filenaming if that differs from the modelname
	# used in querying the file which is the models list variable
	all_models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ]
	modelnames = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4' ]

	modelnames = dict( zip( all_models, modelnames ) )
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )

	os.chdir( output_dir )

	for variable, model, scenario in itertools.product( variables, models, scenarios ):
		# fix the climatology -- precip only
		if variable == 'pr':
			fix_clim = True
		else:
			fix_clim = False
		
		modelname = modelnames[ model ]
		# SETUP BASELINE
		clim_path = os.path.join( base_dir, 'prism', variable )
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
			historical = downscale.Dataset( fn, variable, model, scenario, project=project, units=units, metric=metric, begin=1860, end=2005 )
			future = None # no need for futures here....
		else:
			# get the historical data for anomalies
			historical_fn, = glob.glob( os.path.join( os.path.dirname( fn ).replace( scenario, 'historical' ), '*.nc' ) )
			historical = downscale.Dataset( historical_fn, variable, model, scenario, project=project, units=units, metric=metric, begin=1860, end=2005 )
			future = downscale.Dataset( fn, variable, model, scenario, project=project, units=units, metric=metric, begin=2006, end=2100 )

		# convert from Kelvin to Celcius
		if variable != 'pr':
			if historical:
				historical.ds[ variable ] = historical.ds[ variable ] - 273.15
				historical.ds[ variable ][ 'units' ] = units
			
			if future:
				future.ds[ variable ] = future.ds[ variable ] - 273.15
				future.ds[ variable ][ 'units' ] = units

		if variable == 'pr':
			# convert to mm/month
			if historical:
				timesteps, = historical.ds.time.shape # this assumes time begins in January
				days = [31,28,31,30,31,30,31,31,30,31,30,31] * (timesteps / 12)

				for index, days_in_month in zip(range( len( days ) ), days ):
					historical.ds[ variable ][index, ...] = historical.ds[ variable ][index, ...].data * 86400 * days_in_month

				historical.ds[ variable ][ 'units' ] = units
			
			if future:
				timesteps, = future.ds.time.shape # this assumes time begins in January
				days = [31,28,31,30,31,30,31,31,30,31,30,31] * (timesteps / 12)

				for index, days_in_month in enumerate( days ):
					future.ds[ variable ][index, ...] = future.ds[ variable ][index, ...] * 86400 * days_in_month
					
				future.ds[ variable ][ 'units' ] = units

		# DOWNSCALE
		mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
		clim_begin = '1961'
		clim_end = '1990'

		if variable == 'pr':
			# truncate to whole number
			rounder = np.rint
			downscaling_operation = 'mult'
			aoi_mask = aoi_mask_fn
			# make AOI_Mask input resolution for computing 95th percentiles...
			if aoi_mask_fn is not None:
				aoi_mask = Mask( aoi_mask_fn, historical, 1, 0 )
			else:
				aoi_mask = None
		else:
			# round to 2 decimals
			rounder = partial( np.round, decimals=1 )
			downscaling_operation = 'add'
			aoi_mask = None

		def round_it( x, mask ):
			arr = np.ma.masked_array( data=x, mask=mask )
			return rounder( arr )

		round_data = partial( round_it, mask=( mask==0 ) )

		ar5 = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, future, 
				downscaling_operation=downscaling_operation, mask=mask, mask_value=0, ncpus=32, 
				src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
				post_downscale_function=round_data, varname=variable, modelname=modelname, anom=anom,
				fix_clim=fix_clim, aoi_mask=aoi_mask )

		ar5.downscale( output_dir=output_path )
