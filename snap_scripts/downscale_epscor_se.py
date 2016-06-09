# downscale the prepped cmip5 data downloaded using SYNDA for EPSCoR SE project
# author: Michael Lindgren -- June 09, 2016

import glob, os, rasterio, itertools
import downscale
from downscale import preprocess

# SETUP BASELINE
clim_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/tmax'
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = downscale.Baseline( filelist )

# some setup args
base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prepped_cmip5'
output_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5'
variables = [ 'tasmax', 'tasmin' ]
scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ]

# open a log file to find out where we are messing up
log = open( os.path.join( output_dir, 'log_file_downscale.txt'), 'w' )

for variable, model, scenario in itertools.product( variables, models, scenarios ):
	input_path = os.path.join( base_dir, model, scenario, variable )
	output_path = os.path.join( output_dir, model, scenario, variable )

	if not os.path.exists( output_path ):
		os.makedirs( output_path )

	print( input_path )

	# list files for this set of downscaling -- one per folder
	fn, = glob.glob( os.path.join( input_path, '*.nc' ) )

	if 'historical' in scenario:
		historical = downscale.Dataset( fn, variable, model, scenario, units=None )
		future = None # no need for futures here....
	else:
		# get the historical data for anomalies
		historical_fn, = glob.glob( os.path.join( os.path.dirname( fn ).replace( scenario, 'historical' ), '*.nc' ) )
		historical = downscale.Dataset( historical_fn, variable, model, scenario, units=None )
		future = downscale.Dataset( fn, variable, model, scenario, units=None )

	# DOWNSCALE
	mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
	clim_begin = '1961'
	clim_end = '1990'

	ar5 = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, future, \
			metric='mean', downscaling_operation='add', mask=mask, mask_value=0, ncpus=32, \
			src_crs={'init':'epsg:4326'}, src_nodata=-9999.0, dst_nodata=None,
			post_downscale_function=None )

	ar5.downscale( output_dir=output_path )
