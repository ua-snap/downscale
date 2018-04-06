# RUN THE DOWNSCALER FROM THE COMMAND LINE

if __name__ == '__main__':
	import glob, os, itertools, rasterio
	from downscale import DeltaDownscale, Baseline, Dataset, utils
	import argparse

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-cru", "--cru_ts", action='store', dest='cru_ts', type=str, help="path to the CRU TS3.x NC file to downscale" )
	parser.add_argument( "-o", "--output_path", action='store', dest='output_path', type=str, help="path to the output directory" )
	parser.add_argument( "-clim", "--clim_path", action='store', dest='clim_path', type=str, help="path to the climatology directory -- can only have the climatologies stored in it." )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="cmip5 model name (exact)" )
	parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="cmip5 scenario name (exact)" )
	parser.add_argument( "-p", "--project", action='store', dest='project', type=str, help="cmip5 project name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-ov", "--out_varname", action='store', dest='out_varname', type=str, help="variable name to give to the output files if different from the one stored in the NC file" )
	parser.add_argument( "-met", "--metric", action='store', dest='metric', type=str, help="string name of the metric - mean, max, min, total" )
	parser.add_argument( "-u", "--units", action='store', dest='units', type=str, help="cmip5 units name (exact)" )
	parser.add_argument( "-cb", "--clim_begin", action='store', dest='clim_begin', type=str, help="year of the beginning of the climatology -- ex. '1961'" )
	parser.add_argument( "-ce", "--clim_end", action='store', dest='clim_end', type=str, help="year of the ending of the climatology -- ex. '1990'" )
	parser.add_argument( "-nc", "--ncpus", action='store', dest='ncpus', type=int, help="number of cpus to use in multiprocessing" )
	parser.add_argument( "-do", "--downscaling_operation", action='store', dest='downscaling_operation', type=str, help="downscaling operation to perform.  one of ['add', 'mult']" )
	parser.add_argument( "-i", "--interp", action='store', dest='interp', type=bool, help="boolean of whether to interpolate across space first (True) or not (False) -- used for CRU data" )
	args = parser.parse_args()

	# # unpack the arguments
	cru_ts = args.cru_ts
	output_path = args.output_path
	clim_path = args.clim_path
	model = args.model
	scenario = args.scenario
	project = args.project
	variable = args.variable
	out_varname = args.out_varname
	metric = args.metric
	units = args.units
	clim_begin = args.clim_begin
	clim_end = args.clim_end
	ncpus = args.ncpus
	downscaling_operation = args.downscaling_operation

	# RUN 2.0
	filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
	filelist = [ fn for fn in filelist if not '_14_' in fn ]
	baseline = Baseline( filelist )

	# DOWNSCALE
	mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )

	# FOR CRU WE PASS THE interp=True so we interpolate across space first when creating the Dataset()
	historical = Dataset( cru_ts, variable, model, scenario, project, units, metric=metric, interp=True, method='linear', ncpus=32 )

	ar5 = DeltaDownscale( baseline, clim_begin, clim_end, historical, future=None, \
			downscaling_operation=downscaling_operation, mask=mask, mask_value=0, ncpus=ncpus, \
			src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
			post_downscale_function=None, varname=out_varname, modelname=None ) # -9999

	if not os.path.exists( output_path ):
		os.makedirs( output_path )

	ar5.downscale( output_dir=output_path )
