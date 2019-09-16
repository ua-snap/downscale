# downscale cru data in a CLI way

if __name__ ==	'__main__':
	import glob, os, itertools, rasterio
	from downscale import DeltaDownscale, Baseline, Dataset, utils, Mask
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

	# standard args
	clim_begin = '01-1961'
	clim_end = '12-1990'
	scenario = 'historical'
	project = 'cru'
	anom = True # write out anoms (True) or not (False)
	interp = True # interpolate across space -- Low Res

	# AOI MASK -- HARDWIRE -- GCLL for CRU
	aoi_mask_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/templates/akcan_2km/akcan_aoi_mask_GCLL.shp'

	# RUN 2.0
	filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
	filelist = [ i for i in filelist if '_14_' not in i ] # remove the GD ANNUAL _14_ file.
	baseline = Baseline( filelist )

	# DOWNSCALE
	mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
	
	historical = Dataset( cru_ts, variable, model, scenario, project, units, metric, 
							method='linear', ncpus=32 )

	# post_downscale_function -- rounding
	if variable == 'pr' or variable == 'pre':
		print(variable)
		rounder = np.rint
		downscaling_operation = 'mult'
		find_bounds = True
		fix_clim = True
		# make AOI_Mask at input resolution for computing 95th percentiles...
		if aoi_mask_fn is not None:
			aoi_mask = Mask( aoi_mask_fn, historical, 1, 0 )
		else:
			aoi_mask = None
	else:
		rounder = partial( np.around, decimals=1 )
		downscaling_operation = 'add'
		find_bounds = False
		fix_clim = False
		aoi_mask = None

	def round_it( arr ):
		return rounder( arr )

	# FOR CRU WE PASS THE interp=True so we interpolate across space first when creating the Dataset()
	ar5 = DeltaDownscale( baseline, clim_begin, clim_end, historical, future=None,
				downscaling_operation=downscaling_operation, mask=mask, mask_value=0, ncpus=32,
				src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
				post_downscale_function=round_it, varname=out_varname, modelname=None, 
				anom=anom, interp=interp, find_bounds=find_bounds, fix_clim=fix_clim, aoi_mask=aoi_mask )

	ar5.downscale( output_dir=output_path )
