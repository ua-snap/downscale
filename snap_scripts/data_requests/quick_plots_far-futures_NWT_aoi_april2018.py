# NWT Decadal Precip Trends...
def split_fn( fn ):
	''' split SNAP-formatted files. to Month and Year '''
	dirname, basename = os.path.split( fn )
	basename, ext = os.path.splitext( basename )
	month, year = basename.split( '_' )[ -2: ]
	return { 'month':month, 'year':year, 'fn':fn }

def make_files_df( base_dir, model, scenario, variable, metric='mean' ):
	files = glob.glob( os.path.join( base_dir, model, scenario, variable, '*{}*.tif'.format(metric) ) )
	return pd.DataFrame([ split_fn( fn ) for fn in files ]).sort_values(['year', 'month'])

def global_mean( fn, nodata=-9999 ):
	with rasterio.open( fn ) as rst:
		arr = rst.read( 1 ).copy()
	del rst; rst = None # hard cleanup
	return np.mean(arr[arr != nodata])

if __name__ == '__main__':
	import os, glob
	import matplotlib
	matplotlib.use('agg')
	import rasterio
	import numpy as np
	import pandas as pd
	import multiprocessing as mp
	from matplotlib import pyplot as plt

	# pathing
	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/project_data_delivery/NWT_DELIVERABLES/derived/grids/monthly_decadals'
	output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/project_data_delivery/BrianSieben_NWT_april2018_plots'
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4' ]
	# models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4', 'CRU-TS40' ]
	scenarios = ['rcp45','rcp60','rcp85']
	# scenario = 'historical'
	variables = ['tas','pr']

	for scenario in scenarios:
		for variable in variables:
			out = []
			for model in models:
				print( model )

				# grab filenames to get global area wide means
				df = make_files_df( base_dir, model, scenario, variable )

				# parallel on the viz nodes on ATLAS
				pool = mp.Pool( 64 )
				done = pool.map( global_mean, df.fn.tolist() )
				pool.close()
				pool.join()

				# colname = '{}_{}_{}'.format(model,scenario,variable)
				colname = model # try this instead
				cur_df = pd.DataFrame({
							colname:done}, 
							index=df.apply( lambda x: '{}-{}'.format(x.year, x.month), axis=1))

				years = np.array([ i.split('-')[0] for i in cur_df.index ])
				cur_df = cur_df[ years != '2000s' ]
				year, month = cur_df[~np.isnan(cur_df[colname])].index.tolist()[-1].split('-')
				years = np.array([ i.split('-')[0] for i in cur_df.index ]) # get the years AGAIN
				# futures only here...
				if year == '2100s':
					cur_df = cur_df[ years != '2100s' ]
				elif year == '2300s':
					cur_df = cur_df[ years != '2300s' ]

				out.append( cur_df )

			# concat and write it to disk
			new_df = pd.concat( out, axis=1 )
			# new_df.to_csv( output_filename, sep=',' )

			# now lets plot it.
			variable_lookup = {'pr':'Precipitation','tas':'Temperature'}
			month_lookup = {'01':'January','02':'February','03':'March','04':'April',
							'05':'May','06':'June','07':'July','08':'August',
							'09':'September','10':'October','11':'November','12':'December' }
			label_lookup = {'pr':'millimeters', 'tas':'Celcius'}

			months = ['01','02','03','04','05','06','07','08','09','10','11','12']
			mon_arr = np.array([ i.split('-')[1] for i in new_df.index ])
			for month in months:
				print(month)
				mon_df = new_df[mon_arr == month]
				mon_df.index = [ int(i.split('-')[0][:4]) for i in mon_df.index.tolist() ]

				# title
				title = 'Decadal Mean {} -- {}\nNorthwest Territories\nCMIP5 - {}'.format(variable_lookup[variable], month_lookup[month], scenario)
				# axis title
				# x-axis tick labels
				# standard colors -- lookup table 

				ax = mon_df.plot( kind='line', title=title, figsize=(15,9) )
				ax.set_ylabel( label_lookup[ variable ] )
				ax.set_xlabel( 'Decade' )

				plt.savefig( os.path.join( output_path, 'decadal_mean_monthly_{}_{}_allmodels_month-{}.png'.format(variable,scenario,month)))

				plt.close()

