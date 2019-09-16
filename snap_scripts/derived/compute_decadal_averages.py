# compute decadal averages -- Delta Downscaled SNAP data

def sort_files_dataframe( files ):
	elems = ['variable', 'metric', 'units', 'group', 'model', 'scenario', 'month', 'year']
	# handle CRU naming...  no group...
	if 'CRU-TS40' in files[0]:
		elems = ['variable', 'metric', 'units', 'model', 'scenario', 'month', 'year']
	
	split_fn = [ dict(zip(elems,os.path.basename(fn).split('.tif')[0].split('_') )) for fn in files ]
	_ = [ sfn.update(fn=fn) for sfn, fn in zip(split_fn, files) ]
	df = pd.DataFrame(split_fn)[elems+['fn']]
	return df.sort_values(['year','month']).reset_index(drop=True)
	
def open_raster( fn, band=1 ):
	with rasterio.open( fn ) as rst:
		arr = rst.read( band )
	return arr

def make_decadal( x, output_path ):
	monyear, group = x
	month, decade = monyear
	files = group.fn.tolist()
	arr = np.mean([ open_raster( fn ) for fn in files ], axis=0)
	template = group.iloc[0]

	with rasterio.open( template.fn ) as tmp:
		meta = tmp.meta.copy()
		meta.update(compress='lzw')

	if template['variable'] in ['clt','hurs','tas','tasmax','tasmin',]:
		sig_digits = 1
	elif template['variable'] in ['pr',]:
		sig_digits = 0
	elif template['variable'] in ['vap','rsds',]:
		sig_digits = 4

	out_fn = template.fn.replace( template.year, decade )
	out_fn = os.path.join( output_path, os.path.basename(out_fn) )

	try:
		if not os.path.exists( output_path ):
			_ = os.makedirs( output_path )
	except:
		pass

	with rasterio.open( out_fn, 'w', **meta ) as out:
		arr = np.round(arr, sig_digits)
		out.write( arr, 1 )

	return out_fn


if __name__ == '__main__':
	import os, glob, functools
	import pandas as pd
	import numpy as np
	import rasterio
	import multiprocessing as mp
	import argparse

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='make decadal averages from the SNAP downscaled 2km data' )
	parser.add_argument( "-p", "--path", action='store', dest='path', type=str, help="path to the folders storing the downscaled geotiffs" )
	parser.add_argument( "-o", "--output_path", action='store', dest='output_path', type=str, help="path to the directory to put the decadal averaged output geotiffs" )
	parser.add_argument( "-n", "--ncpus", action='store', dest='ncpus', type=int, help="number of CPU workers to use" )
	args = parser.parse_args()

	# unpack the args
	path = args.path
	output_path = args.output_path
	ncpus = args.ncpus

	# # # # # BEGIN TESTING
	# # this is the full path to the data to perform the decadal on...
	# path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/CRU-TS40/historical/hurs'
	# output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_decadals/CRU-TS40/historical/hurs'
	# ncpus = 32
	# # # # # END TESTING

	ext = '.tif' # hardwire

	# list / sort filenames
	files = glob.glob( os.path.join(path, '*'+ext) )
	files_df = sort_files_dataframe( files )

	# group by decade
	files_df['decade'] = files_df.year.apply( lambda x: x[:3]+'0s' )
	grouped = files_df.groupby(['month','decade'])
	grouped = [ (i,j) for i,j in grouped if len(j) == 10 ]
	
	run = functools.partial( make_decadal, output_path=output_path )
	pool = mp.Pool( ncpus )
	pool.map( run, grouped )
	pool.close()
	pool.join()
