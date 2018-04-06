# downscale the prepped cmip5 data downloaded using SYNDA for EPSCoR SE project
# author: Michael Lindgren -- June 09, 2016

import glob, os, rasterio, itertools
import downscale
from downscale import preprocess

# some setup args
base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prepped_cmip5'
output_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5'
variables = [ 'tasmin', 'tasmax' ]
scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
models = [ 'IPSL-CM5A-LR' ] # , 'MRI-CGCM3''GISS-E2-R', 'GFDL-CM3', 'CCSM4'
# [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ]

# # # 
# modelnames is simply the string name to put in the output filenaming if that differs from the modelname
# used in querying the file which is the models list variable
all_models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ] # temp for distributed run
modelnames = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4' ]
modelnames = dict( zip( all_models, modelnames ) )
# # #

project = 'ar5'
units = 'C'
metric = 'mean'

if not os.path.exists( output_dir ):
	os.makedirs( output_dir )

# open a log file to find out where we are messing up
log = open( os.path.join( output_dir, 'log_file_downscale.txt'), 'w' )

for variable, model, scenario in itertools.product( variables, models, scenarios ):
	modelname = modelnames[ model ]
	# SETUP BASELINE
	if variable == 'tasmin':
		v = 'tmin'
	elif variable == 'tasmax':
		v = 'tmax'
	else:
		NotImplementedError( 'only tasmin and tasmax are currently supported' )

	clim_path = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2', v )
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
		historical = downscale.Dataset( fn, variable, model, scenario, project=project, units=units )
		future = None # no need for futures here....
	else:
		# get the historical data for anomalies
		historical_fn, = glob.glob( os.path.join( os.path.dirname( fn ).replace( scenario, 'historical' ), '*.nc' ) )
		historical = downscale.Dataset( historical_fn, variable, model, scenario, project=project, units=units )
		future = downscale.Dataset( fn, variable, model, scenario, project=project, units=units )

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
			metric=metric, downscaling_operation='add', mask=mask, mask_value=0, ncpus=32, \
			src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
			post_downscale_function=None, modelname=modelname ) # -9999.0

	ar5.downscale( output_dir=output_path )
