# zip up the data for CKAN_Data
def sort_files( files ):
	elems = ['variable', 'metric', 'units', 'group', 'model', 'scenario', 'month', 'year']
	# handle CRU naming...  no group...
	if 'CRU-TS40' in files[0]:
		elems = ['variable', 'metric', 'units', 'model', 'scenario', 'month', 'year']
	
	split_fn = [ dict(zip(elems,os.path.basename(fn).split('.tif')[0].split('_') )) for fn in files ]
	_ = [ sfn.update(fn=fn) for sfn, fn in zip(split_fn, files) ]
	df = pd.DataFrame(split_fn)[elems+['fn']]
	return df.sort_values(['year','month']).fn.tolist()
	
def zip_it( args ):
	dirname, out_fn = args
	return os.system( 'zip -q -rj9 {} {}'.format( out_fn, os.path.join(dirname,'*') ) )

if __name__ == '__main__':
	import os, glob, itertools
	import pandas as pd
	import multiprocessing as mp

	base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_zips_for_CKAN/SNAP'

	variables = ['tas', 'pr', 'clt', 'tasmin', 'tasmax', 'rsds', 'vap', 'hurs',]
	future_models = ['GFDL-CM3', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'NCAR-CCSM4', '5ModelAvg',]
	# future_models = ['5ModelAvg',]
	observed_models = ['CRU-TS40',]

	models = future_models + observed_models

	args = list()
	for model in models:
		# get the right scenarios
		if model == 'CRU-TS40':
			scenarios = ['historical',]
		else:
			scenarios = ['historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85',]

		for scenario, variable in itertools.product(scenarios,variables):
			cur_path = os.path.join( base_path, model, scenario, variable )
			files = glob.glob( os.path.join( cur_path, '*.tif' ) )
			files = sort_files( files ) # sort them
			begin = '_'.join(os.path.basename(files[0]).split('.tif')[0].split('_')[-2:])
			end = '_'.join(os.path.basename(files[-1]).split('.tif')[0].split('_')[-2:])
			out_fn = os.path.join( output_path, '{}_AK_CAN_2km_{}_{}_{}-{}.zip'.format(variable, model, scenario, begin, end ) )
			args = args + [[ os.path.dirname(files[0]), out_fn ]]

	# now zip the data in parallel
	print( 'running...' )
	pool = mp.Pool( 40 )
	out = pool.map( zip_it, args )
	pool.close()
	pool.join()
