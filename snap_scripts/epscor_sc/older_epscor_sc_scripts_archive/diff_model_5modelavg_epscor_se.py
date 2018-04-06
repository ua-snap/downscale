def run( x ):
	return rasterio.open(x).read(1)

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
	years = [ int(fn.split('.')[0].split( split_on )[elem_year]) for fn in files ]
	df = pd.DataFrame( { 'fn':files, 'year':years } )
	df_slice = df[ (df.year >= begin ) & (df.year <= end ) ]
	return df_slice.fn.tolist()


if __name__ == '__main__':
	import numpy as np
	import pandas as pd
	import os, glob, rasterio
	from multiprocessing import Pool

	# MODEL
	input_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5/CCSM4/historical/tasmax'
	
	# list the data
	l = sort_files( glob.glob( os.path.join( input_path, '*.tif' ) ) )
	l = only_years( l, begin=1901, end=2005 )

	# open a pool and turn the list of arrays into an ndarray
	pool = Pool( 31 )
	arr = np.array( pool.map( run, l ) )
	pool.close()
	pool.join()

	# mask it
	arr = np.ma.masked_where( arr <= np.min( arr ), arr )

	# BASELINE -- 5ModelAvg
	input_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5/5ModelAvg/historical/tasmax'
	
	# list the data
	l = sort_files( glob.glob( os.path.join( input_path, '*.tif' ) ) )
	l = only_years( l, begin=1901, end=2005 )

	# open a pool and turn the list of arrays into an ndarray
	pool = Pool( 31 )
	base = np.array( pool.map( run, l ) )
	pool.close()
	pool.join()

	# mask it
	base = np.ma.masked_where( base <= np.min( base ), base )

	# do some stuff with the outputs here
	diff_masked = np.diff( arr - base, axis=0 )

