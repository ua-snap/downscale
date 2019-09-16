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

def diff_mean_domain( baseline, modeled, mask ):
	'''
	baseline = path to 5ModelAvg
	model = path to the single model for same experiment as baseline

	RETURN:
	-------
	mean across space of (mask==1) difference between avg and the model

	'''
	baseline = rasterio.open( baseline ).read( 1 )
	modeled = rasterio.open( modeled ).read( 1 )
	diff = modeled - baseline
	return np.mean( diff[ np.where( mask == 1 ) ] )

def wrap( x ):
	return diff_mean_domain( **x )

if __name__ == '__main__':
	import os, glob, itertools, rasterio
	import xarray as xr
	import numpy as np
	import pandas as pd
	from pathos import multiprocessing as mp

	# setup args
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5'
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_outputs'
	ncpus = 32
	variables = [ 'tasmax' ] # , 'tasmin' ]
	scenarios = [ 'historical' ] # [ 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4', '5ModelAvg' ]
	# decades = [(2010, 2019),(2020, 2029),(2030, 2039),(2040, 2049),(2050, 2059),(2060, 2069),(2070, 2079),(2080, 2089),(2090, 2099)]
	decades = [(1900,1909),(1910, 1919),(1919, 1929),(1930, 1939),(1940, 1949),(1950, 1959),(1960, 1969),(1970, 1979),(1980, 1989),(1990, 1999)]
	template_rst = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5/CCSM4/historical/tasmax/tasmax_mean_C_ar5_CCSM4_historical_1_1901.tif'
	rst = rasterio.open( template_rst )
	subdomain_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/Kenai_StudyArea.shp'

	# create the rasterized version of the input shapefile for spatial query
	subdomains = SubDomains( subdomain_fn, rst, id_field='OBJECTID', name_field='OBJECTID', background_value=0 )
	mask, = subdomains.sub_domains

	# make sure no NoData pixels are in the domain
	mask2 = rst.read_masks(1)
	mask[ mask2 == 0 ] = 0

	all_data = {}
	for variable, model, scenario, decade in itertools.product( variables, models, scenarios, decades ):
		decade_begin, decade_end = decade
		if scenario == 'historical':
			begin = 1900
			end = 2005
		else:
			begin = 2006
			end = 2100

		# mod
		modeled_files = glob.glob( os.path.join( base_path, model, scenario, variable, '*.tif' ) )
		# modeled_files = sort_files( only_years( modeled_files, begin=begin, end=end, split_on='_', elem_year=-1 ) )
		modeled_files = sort_files( only_years( modeled_files, begin=decade_begin, end=decade_end, split_on='_', elem_year=-1 ) )
		
		# groupby month here
		month_grouped = pd.Series( modeled_files ).groupby([ os.path.basename(i).split('_')[-2] for i in modeled_files ])
		month_grouped = { i:j.tolist() for i,j in month_grouped } # make a dict
		
		def f( x ):
			''' apply function for multiprocessing.pool '''
			with rasterio.open( x ) as rst:
				arr = rst.read( 1 )
			return arr

		month_dict = {}
		for month in month_grouped:
			args = month_grouped[ month ]
			# print month
			
			# parallel map
			pool = mp.Pool( ncpus )
			arr = np.array( pool.map( f, args ) )
			pool.close()
			pool.join()
			pool.terminate()
			pool = None
			
			# mask it 
			mask_3d = np.broadcast_to( (mask == 0), arr.shape )
			masked = np.ma.masked_array( arr, mask_3d )
			month_dict[ str(month) ] = { 'stdev':str( np.std( masked ) ),
									'mean':str( np.mean( masked ) ),
									'min':str( np.min( masked ) ),
									'max':str( np.max( masked ) ) }
			
		all_data[ '_'.join([ model, scenario, variable, str(decade_begin), str(decade_end) ]) ] = month_dict
		del arr, mask_3d, masked
		
	# write it out to disk
	if not os.path.exists( output_path ):
		os.makedirs( output_path )
	
	prefix = '_'.join([ variable, 'all_models', 'historical', 'DECADALS', str(begin), str(end) ])

	# this is potentially an issue now with the new dimensionality? -- check at next full run.
	import json
	output_filename = os.path.join( output_path, prefix + '.json' )
	with open( output_filename, 'w' ) as out_json:
		json.dump( all_data, out_json )

	# now some panel-y stuff with the output JSON
	panel = pd.Panel( all_data )
	metrics = ['mean','max','min','stdev']
	for metric in metrics:	
		df = panel[ :, metric, : ].T
		df = df[ [ str(i) for i in range(1,12+1) ] ] # sort the months
		# sort the model combos
		# df = df.reindex_axis([ '_'.join([s,v,m]) for v,m,s in itertools.product(scenarios, variables, models) ], 0)
		# strip variable and underscore
		# df.index = [ ' '.join(i.split('_')[:-1]) for i in df.index ]
		output_filename = os.path.join( output_path, prefix + '_' + metric +'.csv' )
		df.to_csv( output_filename, sep=',' )


	# [ '_'.join([v,m,s]) for v,m,s in variables, models, scenarios ]

	# df = pd.DataFrame( all_means )
	# df.to_csv( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived/epscor_se_krw_mean_diffs_model_vs_5ModelAvg_monthly.csv', sep=',' )
