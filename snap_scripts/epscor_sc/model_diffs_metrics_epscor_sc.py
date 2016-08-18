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

def diff( x, y, *args, **kwargs ):
	'''
	difference between 2 np.arrays representing 
	2-D rasters in the format : GeoTiff

	ARGUMENTS:
	----------
	x = [str] path to the baseline GeoTiff
	y = [str] path to the modeled GeoTiff

	RETURNS:
	-------
	difference of the 2 arrays as a 2D numpy array
	( y - x )
	
	'''
	import rasterio
	baseline = rasterio.open( x ).read( 1 )
	modeled = rasterio.open( y ).read( 1 )
	return ( modeled - baseline  )

def wrap_diff( x ):
	''' simple wrapper for multiprocessing '''
	return diff( *x )

def get_metrics( base_path, variable, model, scenario, decade, mask, domain_name=None, ncpus=32 ):
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
		modeled = month_grouped[ month ]
		# change the model name to the baseline model in the series for comparison
		baseline = [ fn.replace( model, '5ModelAvg' ) for fn in modeled ]
		args = zip( baseline, modeled )

		# get diffs in parallel
		pool = mp.Pool( ncpus )
		arr = np.array( pool.map( wrap_diff, args ) )
		pool.close()
		pool.join()
		pool.terminate()
		pool = None

		# this derives a mean from 3D (time, x, y) to 2D (x, y)
		mean_arr = np.mean( arr, axis=0 )
		arr = None
		masked = np.ma.masked_array( mean_arr, mask == 0 )

		# calculate metrics across the 2D space
		month_dict[ str(month) ] = { 'stdev':str( np.std( masked ) ),
									'mean':str( np.mean( masked ) ),
									'min':str( np.min( masked ) ),
									'max':str( np.max( masked ) ) }

		# domain_name
		if domain_name == None:
			domain_name, = str( np.unique( mask > 0 ) )
	
	return { '_'.join([ model, scenario, variable, domain_name, str(decade_begin), str(decade_end) ]) : month_dict }

if __name__ == '__main__':
	import os, glob, itertools, rasterio, json
	from copy import deepcopy
	import xarray as xr
	import numpy as np
	import pandas as pd
	from pathos import multiprocessing as mp

	# setup args
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/EPSCOR_SC_DELIVERY_AUG2016/derived/grids'
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/EPSCOR_SC_DELIVERY_AUG2016/derived/tabular'
	ncpus = 32
	project = 'cmip5'
	variables = [ 'tasmin', 'tasmax', 'tas', 'pr' ]
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4' ]
	scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	template_rst = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/EPSCOR_SC_DELIVERY_AUG2016/downscaled/NCAR-CCSM4/historical/tasmax/tasmax_mean_C_ar5_NCAR-CCSM4_historical_01_1901.tif'
	rst = rasterio.open( template_rst )
	# subdomain_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/Kenai_StudyArea.shp'
	subdomain_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/SCTC_watersheds.shp'
	
	# create the rasterized version of the input shapefile for spatial query
	# subdomains = SubDomains( subdomain_fn, rst, id_field='OBJECTID', name_field='OBJECTID', background_value=0 )
	subdomains = SubDomains( subdomain_fn, rst, id_field='OBJECTID', name_field='HU_12_Name', background_value=0 )
	masks = subdomains.sub_domains

	# make sure no NoData pixels are in the domain
	nodata_mask = rst.read_masks( 1 ) # mask where zero
	for count, mask in enumerate( masks ):
		mask[ nodata_mask == 0 ] = 0
		masks[ count ] = mask

	for variable in variables:
		print( variable )
		all_data = {}
		for model, scenario in itertools.product( models, scenarios ):
			if scenario == 'historical':
				decades = [(1900,1909),(1910, 1919),(1920, 1929),(1930, 1939),(1940, 1949),\
							(1950, 1959),(1960, 1969),(1970, 1979),(1980, 1989),(1990, 1999),(2000,2005)]
			else:
				decades = [(2006,2009),(2010, 2019),(2020, 2029),(2030, 2039),(2040, 2049),(2050, 2059),\
							(2060, 2069),(2070, 2079),(2080, 2089),(2090, 2099)]

			for decade in decades:
				if scenario == 'historical':
					begin = 1900
					end = 2005
				else:
					begin = 2006
					end = 2100
			
				print( 'running: {} {} {} {}'.format( model, variable, scenario, decade ) )

				for mask in masks:
					domain_num, = np.unique(mask[mask > 0])
					domain_name = subdomains.names_dict[ domain_num ].replace( ' ', '' )
					# run it

					all_data.update( get_metrics( base_path, variable, model, scenario, decade, mask, domain_name, ncpus ) )
		
		# write it out to disk
		if not os.path.exists( output_path ):
			os.makedirs( output_path )
		
		prefix = '_'.join([ variable, project, 'decadals', '5ModelAvg_diff','summaries', str(1900), str(2100) ])

		# its LONG FORMAT output with all datas in rows for a single var/metric
		output_filename = os.path.join( output_path, prefix + '.json' )
		with open( output_filename, 'w' ) as out_json:
			json.dump( all_data, out_json )

		# now some panel-y stuff with the output JSON
		panel = pd.Panel( deepcopy( all_data ) ).copy()
		metrics = ['mean','max','min','stdev']
		for metric in metrics:
			df = panel[ :, metric, : ].T
			# sort the months
			df = df[ [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ] ]
			output_filename = os.path.join( output_path, prefix + '_' + metric +'.csv' )
			if variable == 'pr':
				df = df.astype( np.float ).round( 0 ).astype( np.int )
				# output to csv -- int so no rounding needed.
				df.to_csv( output_filename, sep=',')
			else:
				df = df.astype( np.float32 )
				# round tas* to single decimal place
				df = df.copy().round( decimals=1 )
				# output ensuring single decimal place as string
				df.apply( lambda x: x.apply( lambda x: '%2.1f' % x) ).to_csv( output_filename, sep=',', float_format='%11.6f')
