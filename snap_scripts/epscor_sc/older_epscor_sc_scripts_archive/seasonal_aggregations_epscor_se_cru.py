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
def coordinates( fn=None, meta=None, numpy_array=None, input_crs=None, to_latlong=False ):
	'''
	take a raster file as input and return the centroid coords for each 
	of the grid cells as a pair of numpy 2d arrays (longitude, latitude)
	'''
	import rasterio
	import numpy as np
	from affine import Affine
	from pyproj import Proj, transform

	if fn:
		# Read raster
		with rasterio.open( fn ) as r:
			T0 = r.affine  # upper-left pixel corner affine transform
			p1 = Proj( r.crs )
			A = r.read( 1 )  # pixel values

	elif (meta is not None) & (numpy_array is not None):
		A = numpy_array
		if input_crs != None:
			p1 = Proj( input_crs )
			T0 = meta[ 'affine' ]
		else:
			p1 = None
			T0 = meta[ 'affine' ]
	else:
		BaseException( 'check inputs' )

	# All rows and columns
	cols, rows = np.meshgrid(np.arange(A.shape[1]), np.arange(A.shape[0]))
	# Get affine transform for pixel centres
	T1 = T0 * Affine.translation( 0.5, 0.5 )
	# Function to convert pixel row/column index (from 0) to easting/northing at centre
	rc2en = lambda r, c: ( c, r ) * T1
	# All eastings and northings (there is probably a faster way to do this)
	eastings, northings = np.vectorize(rc2en, otypes=[np.float, np.float])(rows, cols)
 	
	if to_latlong == False:
		return eastings, northings
	elif (to_latlong == True) & (input_crs != None):
		# Project all longitudes, latitudes
		longs, lats = transform(p1, p1.to_latlong(), eastings, northings)
		return longs, lats
	else:
		BaseException( 'cant reproject to latlong without an input_crs' )

def get_month_seaon( fn ):
	# seasons
	seasonal_lookup = { 1:'DJF', 2:'DJF', 3:'MAM', 4:'MAM', 5:'MAM', \
						6:'JJA', 7:'JJA', 8:'JJA',\
						 9:'SON', 10:'SON', 11:'SON', 12:'DJF' }

	fn = os.path.basename( fn )
	month, year = fn.replace( '.tif', '' ).split( '_' )[-2:]
	return seasonal_lookup[ int(month) ]

def get_year( fn ):
	fn = os.path.basename( fn )
	month, year = fn.replace( '.tif', '' ).split( '_' )[-2:]
	return year

# def f( x ):
# 	import rasterio
# 	season = get_month_seaon( x[0] )
# 	meta = rasterio.open( x[0] ).meta
# 	meta.update( compress='lzw' )
# 	mean = np.mean([ rasterio.open(i).read(1) for i in x ], axis=0 )
# 	dirname, basename = os.path.split( x[0] )
# 	bsplit = basename.split( '.' )[0].split( '_' )
# 	bsplit[-2] = season
# 	output_path = dirname.replace( '/downscaled_', 'SEASONALS/downscaled_' )
# 	try:
# 		if not os.path.exists( output_path ):
# 			os.makedirs( output_path )
# 	except:
# 		pass
# 	output_filename = os.path.join( output_path, '_'.join( bsplit ) + '.tif' )
# 	with rasterio.open( output_filename, 'w', **meta ) as out:
# 		out.write( mean, 1 )

def to_gtiff( x, meta, output_path, prefix, *args, **kwargs ):
	import rasterio
	decade_season, arr = x
	decade, season = decade_season
	
	try:
		if not os.path.exists( output_path ):
			os.makedirs( output_path )
	except:
		pass

	output_filename = os.path.join( output_path, '_'.join( [prefix, season, decade + '.tif'] ) )
	with rasterio.open( output_filename, 'w', **meta ) as out:
		arr[ np.isinf( arr ) ] = meta[ 'nodata' ]
		out.write( arr, 1 )
	return output_filename

def wrap( x ):
	''' 
	multiprocessing wrapper for clean 
	argument handling without lambda 
	'''
	return to_gtiff( *x )

def main( base_path, output_path, model, scenario, variable, begin, end, ncpus ):
	'''
	function to calculate and output mean seasonal monthly data across decades
	
	ARGUMENTS:
	----------
	base_path = [  ]  
	output_path = [  ]  
	model = [  ]  
	scenario = [  ]  
	variable = [  ]  
	begin = [  ]  
	end = [  ]  
	ncpus = [  ]  

	RETURNS
	-------
	output_directory of newly produced GeoTiffs if successful. else, error.

	'''
	# modeled data
	files = glob.glob( os.path.join( base_path, model, scenario, variable, '*.tif' ) )
	files = sort_files( only_years( files, begin=begin, end=end, split_on='_', elem_year=-1 ) )
	decade_grouper = [ get_year(fn)[:3]+'0' for fn in files ]
	decades = pd.Series( files ).groupby( decade_grouper )

	# THIS CURRENTLY WORKS AND GIVES SEASONAL AVERAGES FOR EVERY DECADE IN THE SERIES!
	seasonal_decade_means = decades.apply( lambda x: x.groupby( [ get_month_seaon( i ) for i in x ] ) \
				.apply( lambda y: np.mean([ rasterio.open(i).read(1) for i in y ], axis=0 ) ) )

	# set up metadata
	meta = rasterio.open( files[0] ).meta
	if 'transform' in meta.keys():
		meta.pop( 'transform' )
	meta.update( compress='lzw' )

	# output_prefix
	prefix = '_'.join( os.path.basename( files[0] ).split( '_' )[:-2] )

	# now we should output these to disk...
	output_path_join = os.path.join( output_path, model, scenario, variable )
	
	pool = mp.Pool( ncpus )
	out = pool.map( wrap, [( i, meta, output_path, prefix ) for i in seasonal_decade_means.iteritems()] )
	pool.close()
	pool.join()
	pool.terminate()
	pool = None
	return output_path_join

if __name__ == '__main__':
	import os, glob, itertools, rasterio
	import xarray as xr
	import pandas as pd
	import numpy as np
	from pathos import multiprocessing as mp

	# some setup
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_clipped'
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_outputs/seasonal_decadals'
	models = ['ts323']# [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4', '5ModelAvg' ]
	scenarios = [ 'historical' ] #, 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	variables = [ 'tasmin', 'tasmax', 'pr', 'tas' ]
	ncpus = 32

	# run all combinations
	for model, scenario, variable in itertools.product(models, scenarios, variables):
		out_path = os.path.join( output_path, model, scenario, variable )
		if scenario == 'historical':
			begin = 1900
			end = 2005
		else:
			begin = 2006
			end = 2100
		print( 'running: {} {} {}'.format( model, scenario, variable ) )
		new_dirs = main( base_path, out_path, model, scenario, variable, begin, end, ncpus )

