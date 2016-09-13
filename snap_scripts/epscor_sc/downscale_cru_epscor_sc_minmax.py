# downscale cru data in a CLI way
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
	begin = 1901
	end = 2014

	# standard args
	clim_begin = '01-1961'
	clim_end = '12-1990'
	scenario = 'historical'
	project = 'cru'
	anom = False # write out anoms (True) or not (False)
	# output_dir = os.path.join( base_dir, 'downscaled_minmax' )
	# output_path = os.path.join( output_dir, model, scenario, variable )

	# RUN 2.0
	# filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
	# # filelist = [ i for i in filelist if '_14_' not in i ] # remove the GD ANNUAL _14_ file.
	# baseline = Baseline( filelist )

	# clim_path = os.path.join( base_dir, 'downscaled', modelname, scenario, mean_variable_out )
	filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
	# sort these files
	filelist = only_years( sort_files( filelist ), begin=begin, end=end )
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

	# FOR CRU WE PASS THE interp=True so we interpolate across space first when creating the Dataset()
	historical = Dataset( cru_ts, variable, model, scenario, project, units, metric, 
							interp=True, method='linear', ncpus=32 )

	mean_fn, = cru_ts.replace( variable, mean_variable )
	mean_ds = downscale.Dataset( mean_fn, mean_variable_cru, model, scenario, project=project, units=units, metric=metric, begin=begin, end=end )

	ar5 = DeltaDownscaleMinMax( baseline, clim_begin, clim_end, historical, future=None,
				downscaling_operation=downscaling_operation, mask=mask, mask_value=0, ncpus=32,
				src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
				post_downscale_function=round_it, varname=out_varname, modelname=None, anom=False, 
				mean_ds=mean_ds, mean_variable=mean_variable_cru )

	if not os.path.exists( output_path ):
		os.makedirs( output_path )

	ar5.downscale( output_dir=output_path )
