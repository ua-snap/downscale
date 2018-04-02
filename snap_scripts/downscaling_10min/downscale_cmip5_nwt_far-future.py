# downscale the prepped cmip5 data downloaded using SYNDA @ 10min
# author: Michael Lindgren -- 2017 TAS/PR DATA RUN...

if __name__ == '__main__':
	import glob, os, rasterio, itertools
	from functools import partial
	from downscale import preprocess, Baseline, Mask, utils, DeltaDownscaleFF, DatasetFF
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
	parser.add_argument( "-by", "--begin", action='store', dest='begin', nargs='?', const=0, type=int, help="desired year to begin downscaling if diff than whole dset" )
	parser.add_argument( "-ey", "--end", action='store', dest='end', nargs='?', const=0, type=int, help="desired year to end downscaling if diff than whole dset" )

	args = parser.parse_args()

	# unpack the args
	variable = args.variable
	scenario = args.scenario
	model = args.model
	units = args.units
	metric = args.metric
	base_dir = args.base_dir
	begin = args.begin
	end = args.end

	if begin == 0:
		begin = None
	
	if end == 0:
		end = None

	# AOI MASK -- HARDWIRED -- PCLL for CMIP5
	aoi_mask_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/masks/pcll_template_10min_extent_with_nwt.shp'
	project = 'ar5'
	
	# # # # FOR TESTING # # # 
	# base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
	# # variable = 'tas'
	# variable = 'pr'
	# scenario = 'rcp85'
	# model = 'NCAR-CCSM4'
	# # units = 'C'
	# units = 'mm'
	# # metric = 'mean'
	# metric = 'total'
	# begin = None
	# end = None
	# # END TESTING

	# some setup args
	base_path = os.path.join( base_dir,'cmip5_nwt','cmip5_raw_ncrcat' ) # NEW PATH FOR THE NCRCAT FILES... WORKS WELL.
	output_dir = os.path.join( base_dir,'downscaled_10min_nwt' )
	variables = [ variable ]
	scenarios = [ scenario ]
	models = [ model ]
	anom = True # write out anoms (True) or not (False)
	interp = False # interpolate across space -- Low Res
	find_bounds = False

	# modelnames is simply the string name to put in the output filenaming if that differs from the modelname
	# used in querying the file which is the all_models list variable
	all_models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4', 'CCSM4' ]
	modelnames = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4', 'NCAR-CCSM4' ]

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
		variable_lookup_cru = {'pr':'pre', 'tas':'tmp'} # only ready for tas / pr currently
		clim_path = os.path.join( base_dir, 'cru', 'akcan_10min_extent','cru_cl20', variable_lookup_cru[ variable ] )
		filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
		baseline = Baseline( filelist )
		
		input_path = os.path.join( base_path, model, scenario, variable )
		output_path = os.path.join( output_dir, model, scenario, variable )

		if not os.path.exists( output_path ):
			os.makedirs( output_path )

		print( input_path )

		# list files for this set of downscaling
		fn = sorted( glob.glob( os.path.join( input_path, '*.nc' ) ) )

		if 'historical' in scenario:
			historical = DatasetFF( fn, variable, model, scenario, project=project, units=units, metric=metric, begin=begin, end=end )
			future = None # no need for futures here....
		else:
			# get the historical data for anomalies
			historical_fn = sorted( glob.glob( os.path.join( os.path.dirname( fn[0] ).replace( scenario, 'historical' ), '*.nc' ) ) )
			# [ NOTE ]: begin/end=None bc want all hist avail
			historical = DatasetFF( historical_fn, variable, model, scenario, project=project, units=units, metric=metric, begin=None, end=None ) 
			future = DatasetFF( fn, variable, model, scenario, project=project, units=units, metric=metric, begin=begin, end=end )

		# convert from Kelvin to Celcius
		if variable != 'pr':
			if historical:
				historical.dat = historical.dat - 273.15
				historical.units = units
			
			if future:
				future.dat = future.dat - 273.15
				future.units = units

		if variable == 'pr':
			# convert to mm/month
			if historical:
				timesteps = historical.dat.shape[0] # this assumes time begins in January
				days = [31,28,31,30,31,30,31,31,30,31,30,31] * int(timesteps / 12)

				for index, days_in_month in zip(range( len( days ) ), days ):
					historical.dat[index, ...] = historical.dat[index, ...] * 86400 * days_in_month

				historical.units = units
			
			if future:
				timesteps = future.dat.shape[0] # this assumes time begins in January
				days = [31,28,31,30,31,30,31,31,30,31,30,31] * int(timesteps / 12)

				for index, days_in_month in enumerate( days ):
					future.dat[index, ...] = future.dat[index, ...] * 86400 * days_in_month
					
				future.units = units

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

		ar5 = DeltaDownscaleFF( baseline, clim_begin, clim_end, historical, future, 
				downscaling_operation=downscaling_operation, mask=mask, mask_value=0, ncpus=32, 
				src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
				post_downscale_function=round_data, varname=variable, modelname=modelname, anom=anom,
				fix_clim=fix_clim, aoi_mask=aoi_mask )

		ar5.downscale( output_dir=output_path )
