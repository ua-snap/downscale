# maybe read in the baseline
# then loop through reads of all models...
# perform the diff
# then groupby month and compute means / stdev

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

class SubDomains( object ):
	'''
	rasterize subdomains shapefile to ALFRESCO AOI of output set
	'''
	def __init__( self, subdomains_fn, rasterio_raster, id_field, name_field, background_value=0, *args, **kwargs ):
		'''
		initializer for the SubDomains object
		The real magic here is that it will use a generator to loop through the 
		unique ID's in the sub_domains raster map generated.
		'''
		import numpy as np
		self.subdomains_fn = subdomains_fn
		self.rasterio_raster = rasterio_raster
		self.id_field = id_field
		self.name_field = name_field
		self.background_value = background_value
		self._rasterize_subdomains( )
		self._get_subdomains_dict( )

	def _rasterize_subdomains( self ):
		'''
		rasterize a subdomains shapefile to the extent and resolution of 
		a template raster file. The two must be in the same reference system 
		or there will be potential issues. 
		returns:
			numpy.ndarray with the shape of the input raster and the shapefile
			polygons burned in with the values of the id_field of the shapefile
		gotchas:
			currently the only supported data type is uint8 and all float values will be
			coerced to integer for this purpose.  Another issue is that if there is a value
			greater than 255, there could be some error-type issues.  This is something that 
			the user needs to know for the time-being and will be fixed in subsequent versions
			of rasterio.  Then I can add the needed changes here.
		'''
		import geopandas as gpd
		import numpy as np

		gdf = gpd.read_file( self.subdomains_fn )
		id_groups = gdf.groupby( self.id_field ) # iterator of tuples (id, gdf slice)

		out_shape = self.rasterio_raster.height, self.rasterio_raster.width
		out_transform = self.rasterio_raster.affine

		arr_list = [ self._rasterize_id( df, value, out_shape, out_transform, background_value=self.background_value ) for value, df in id_groups ]
		self.sub_domains = arr_list
	@staticmethod
	def _rasterize_id( df, value, out_shape, out_transform, background_value=0 ):
		from rasterio.features import rasterize
		geom = df.geometry
		out = rasterize( ( ( g, value ) for g in geom ),
							out_shape=out_shape,
							transform=out_transform,
							fill=background_value )
		return out
	def _get_subdomains_dict( self ):
		import geopandas as gpd
		gdf = gpd.read_file( self.subdomains_fn )
		self.names_dict = dict( zip( gdf[self.id_field], gdf[self.name_field] ) )

def read_raster( x ):
	''' apply function for multiprocessing.pool 
		helps with clean i/o '''
	with rasterio.open( x ) as rst:
		arr = rst.read( 1 )
	return arr

def make_decadals( base_path, output_path, variable, model, scenario, decade, ncpus, agg_metric ):
	'''
	main function to return monthly summary stats for the group
	as a `dict`
	'''
	decade_begin, decade_end = decade
	modeled_files = glob.glob( os.path.join( base_path, model, scenario, variable, '*.tif' ) )
	modeled_files = sort_files( only_years( modeled_files, begin=decade_begin, end=decade_end, split_on='_', elem_year=-1 ) )
	
	# groupby month here
	month_grouped = pd.Series( modeled_files ).groupby([ os.path.basename(i).split('_')[-2] for i in modeled_files ])
	month_grouped = { i:j.tolist() for i,j in month_grouped } # make a dict
	
	month_dict = {}
	for month in month_grouped:
		pool = mp.Pool( ncpus )
		arr = np.array( pool.map( read_raster, month_grouped[ month ] ) )
		pool.close()
		pool.join()
		pool.terminate()
		pool = None

		template = month_grouped[ month ][0]
		var, metric, units, project, model, scenario, month, year = \
			os.path.basename( template ).split( '.' )[0].split( '_' )
		rst = rasterio.open( template )
		meta = rst.meta
		mask = rst.read_masks( 1 )

		if 'transform' in meta:
			meta.pop( 'transform' )
		meta.update( compress='lzw' )

		metric_switch = { 'mean':np.mean, 'total':np.sum, 'min':np.min, 'max':np.max }
		variable, metric, units, project, model, scenario = os.path.basename( fn ).split( '.' )[0].split( '_' )[:-2]
		arr = metric_switch[ metric_agg ]( [ read_raster( i ) for i in files ], axis=0 )

		# this derives a mean from 3D (time, x, y) to 2D (x, y)
		mean_arr = np.mean( arr, axis=0 )
		decade_out = str(decade_begin)[:3] + '0s'
		output_filename = os.path.join( output_path, model, scenario, variable, '_'.join([var, metric_agg, units, project, model, scenario, month, decade_out]) + '.tif' )

		# make sure the dirname exists
		dirname = os.path.dirname( output_filename )
		try:
			if not os.path.exists( dirname ):
				os.makedirs( dirname )
		except:
			pass
		
		with rasterio.open( output_filename, 'w', **meta ) as out:
			mean_arr[ mask == 0 ] = meta[ 'nodata' ]
			out.write( mean_arr, 1 )

	return output_filename


if __name__ == '__main__':
	import os, glob, itertools, rasterio, json
	from copy import deepcopy
	import xarray as xr
	import numpy as np
	import pandas as pd
	from pathos import multiprocessing as mp
	import argparse

	'''
	this tool assumes that the data are stored in a directory structure as follows:
	
	base_path
		model
			scenario
				variable
					FILES
	'''

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-b", "--base_path", action='store', dest='base_path', type=str, help="path to the directory where the downscaled modeled data are stored" )
	parser.add_argument( "-o", "--output_path", action='store', dest='output_path', type=str, help="path to the output directory" )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="cmip5 model name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-am", "--agg_metric", action='store', dest='agg_metric', type=str, help="string name of the metric to compute the decadal summary - mean, max, min, total" )
	parser.add_argument( "-nc", "--ncpus", action='store', dest='ncpus', type=int, help="number of cpus to use in multiprocessing" )	
	args = parser.parse_args()
	
	# unpack for cleaner var access:
	base_path = args.base_path
	output_path = args.output_path
	model = args.model
	variable = args.variable
	ncpus = args.ncpus
	agg_metric = args.agg_metric

	# switches to deal with different date groups.  Hardwired to CMIP5 and CRU TS323 currently.
	cmip_switch = { 'historical':(1900,2005), 'rcp26':(2005,2100), 'rcp45':(2005,2100), 'rcp60':(2005,2100), 'rcp85':(2006,2100) }
	cru_switch = { 'historical':(1901,2014) }
	project_switch = { 'cmip5':cmip_switch, 'cru':cru_switch }

	# decades setup
	begin, end = project_switch[ project ][ scenario ]
	decades = list( sorted( set( [ (int(str(i)[:3]+'0'), int(str(i)[:3]+'9')) for i in range( begin, end ) ] ) ) )

	for decade in decades:
		print( 'running: {} {} {} {}'.format( model, variable, scenario, decade ) )
		_ = make_decadals( base_path, output_path, variable, model, scenario, decade, ncpus, agg_metric )




# FOR TESTING THEN REMOVE
	# # setup args
	# base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_clipped'
	# output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_outputs/annual_decadals'
	# ncpus = 32
	# project = 'cru' # 'cmip5'
	# variables = [ 'tasmin', 'tasmax', 'tas', 'pr' ]
	# models = [ 'ts323' ] # [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4', '5ModelAvg' ]
	# scenarios = [ 'historical' ] # [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ] # 
	# agg_metric = 'mean'




# # THE OLDER WAY TO DO IT:
# 	for variable in variables:
# 		for model, scenario in itertools.product( models, scenarios ):
			
# 			begin, end = project_switch[ project ][ scenario ]
# 			decades = list( sorted( set( [ (int(str(i)[:3]+'0'), int(str(i)[:3]+'9')) for i in range( begin, end ) ] ) ) )

# 			# if scenario == 'historical':
# 			# 	# decades = [(1900,1909),(1910, 1919),(1920, 1929),(1930, 1939),(1940, 1949),\
# 			# 	# 			(1950, 1959),(1960, 1969),(1970, 1979),(1980, 1989),(1990, 1999),(2000,2005)]
# 			# 	decades = [(1901,1909),(1910, 1919),(1920, 1929),(1930, 1939),(1940, 1949),\
# 			# 				(1950, 1959),(1960, 1969),(1970, 1979),(1980, 1989),(1990, 1999),(2000,2009),(2010, 2014)]
# 			# else:
# 			# 	decades = [(2006,2009),(2010, 2019),(2020, 2029),(2030, 2039),(2040, 2049),(2050, 2059),\
# 			# 				(2060, 2069),(2070, 2079),(2080, 2089),(2090, 2099)]

# 			for decade in decades:
# 				# if scenario == 'historical':
# 				# 	if project == 'cru':
# 				# 		begin = 1901
# 				# 		end = 2014
# 				# else:
# 				# 	begin = 2006
# 				# 	end = 2100

# 				print( 'running: {} {} {} {}'.format( model, variable, scenario, decade ) )
				
# 				final = make_decadals( base_path, output_path, variable, model, scenario, decade, ncpus=32 )

