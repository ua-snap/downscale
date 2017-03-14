# downscale cru data in a CLI way
# downscale the prepped cmip5 data downloaded using SYNDA for EPSCoR SC project
# author: Michael Lindgren -- June 09, 2016

if __name__ ==	'__main__':
	import glob, os, itertools, rasterio
	import downscale
	from downscale import DeltaDownscaleMinMax, Baseline, Dataset, utils
	from functools import partial
	import numpy as np
	import argparse

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-ts", "--ts", action='store', dest='cru_ts', type=str, help="path to the cru file to use in downscaling (.nc)" )
	parser.add_argument( "-cl", "--clim_path", action='store', dest='clim_path', type=str, help="path to the directory where the 12 geotiff climatology files are stored" )
	parser.add_argument( "-o", "--output_path", action='store', dest='output_path', type=str, help="path to the output directory" )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="cmip5 model name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-mvc", "--mean_variable_cru", action='store', dest='mean_variable_cru', type=str, help="cru mean variable name (exact)" )
	parser.add_argument( "-mvo", "--mean_variable_out", action='store', dest='mean_variable_out', type=str, help="cru mean variable output name (exact)" )
	parser.add_argument( "-u", "--units", action='store', dest='units', type=str, help="string name of the units data are in" )
	parser.add_argument( "-met", "--metric", action='store', dest='metric', type=str, help="string name of the metric data are in" )	
	parser.add_argument( "-nc", "--ncpus", action='store', dest='ncpus', type=int, help="number of cpus to use in multiprocessing" )
	parser.add_argument( "-ov", "--out_varname", action='store', dest='out_varname', type=str, help="string name of output name to use instead of variable in file" )
	parser.add_argument( "-b", "--begin", action='store', dest='begin', type=int, help="4 digit beginning year (i.e. 1901)" )
	parser.add_argument( "-e", "--end", action='store', dest='end', type=int, help="4 digit ending year (i.e. 2015)" )
	args = parser.parse_args()

	# unpack args
	cru_ts = args.cru_ts
	clim_path = args.clim_path
	output_path = args.output_path
	model = args.model
	variable = args.variable
	units = args.units
	metric = args.metric
	ncpus = args.ncpus
	out_varname = args.out_varname
	mean_variable_cru = args.mean_variable_cru
	mean_variable_out = args.mean_variable_out
	begin = args.begin
	end = args.end
	# begin = 1901
	# end = 2014

	# standard args
	clim_begin = '01-1961'
	clim_end = '12-1990'
	scenario = 'historical'
	project = 'cru'
	anom = True # write out anoms (True) or not (False)
	interp = True

	# clim_path = os.path.join( base_dir, 'downscaled', modelname, scenario, mean_variable_out )
	filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
	# sort these files
	# filelist = utils.only_years( utils.sort_files( filelist ), begin=begin, end=end )
	filelist = utils.sort_files( filelist ) # [CHECK] this changed to remove begin/end vars...  not needed...
	baseline = downscale.Baseline( filelist )

	# DOWNSCALE
	mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )

	# make round/trunc function for post_downscale_function
	if variable == 'pr' or variable == 'pre':
		rounder = np.rint
		downscaling_operation = 'mult'
	else:
		rounder = partial( np.around, decimals=1 )
		downscaling_operation = 'add'

	def round_it( arr ):
		return rounder( arr )

	historical = Dataset( cru_ts, variable, model, scenario, project, units, metric, 
							method='linear', ncpus=32 )

	mean_fn = cru_ts.replace( variable, mean_variable_cru )
	mean_ds = downscale.Dataset( mean_fn, mean_variable_cru, model, scenario, project=project, units=units, metric=metric, begin=begin, end=end )

	# FOR CRU WE PASS THE interp=True so we interpolate across space first when creating the Dataset()
	ar5 = DeltaDownscaleMinMax( baseline=baseline, clim_begin=clim_begin, clim_end=clim_end, historical=historical, future=None,
				downscaling_operation=downscaling_operation, mask=mask, mask_value=0, ncpus=32,
				src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
				post_downscale_function=round_it, varname=out_varname, modelname=None, anom=anom, 
					mean_ds=mean_ds, mean_variable=mean_variable_cru, interp=interp )

	if not os.path.exists( output_path ):
		os.makedirs( output_path )

	ar5.downscale( output_dir=output_path )
