# make 5 model avg of the data in 
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

def generate( files ):
	rst = rasterio.open( files[0] )
	meta = rst.meta
	meta.update( compress='lzw' )
	mask = rst.read_masks( 1 )

	if 'transform' in meta.keys():
		meta.pop( 'transform' )

	arr_group = np.array([ rasterio.open( fn ).read( 1 ) for fn in files ])
	group_mean = np.mean( arr_group, axis=0 )
	group_mean[ mask == 0 ] = -3.39999995e+38 # add back the oob

	fn, = [ i for i in files if 'CCSM4' in i ]
	output_filename = fn.replace( 'CCSM4', '5ModelAvg' )
	new_path = os.path.dirname( output_filename )
	try:
		# new directory -- if needed:
		if not os.path.exists( new_path ):
			os.makedirs( new_path )
	except:
		pass

	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( group_mean, 1 )
	return output_filename

if __name__ == '__main__':
	import os, sys, glob, rasterio, shutil
	import pandas as pd
	import numpy as np
	from pathos import multiprocessing as mp

	base_path = '/Data/Base_Data/Climate/AK_CAN_2km/historical/AR5_CMIP5_models'
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ]
	variables = [ 'pr', 'tas' ]

	for variable in variables:
		files_out = {}
		for model in models:
			out = [ [ os.path.join( root, i ) for i in files ] for root, subs, files in os.walk( os.path.join( base_path, model ) ) \
							if 'derived' not in root \
							if len([ fn for fn in files if fn.endswith('.tif') ]) > 0 \
							if variable in root ]
			out = [ j for i in out for j in i ]
			begin = 1900
			end = 2005
			out_sorted = sort_files( only_years( out, begin=begin, end=end, split_on='_', elem_year=-1 ) )
			files_out.update({ model:out_sorted })

		# put the files in a pd.DataFrame obj
		df = pd.DataFrame( files_out )
		args_list = [ dat.tolist() for row, dat in df.iterrows() ]

	 	pool = mp.Pool( 32 )
	 	out = pool.map( generate, args_list )
	 	pool.close()
	 	pool.join()
	 	pool.terminate()
	 	pool = None

