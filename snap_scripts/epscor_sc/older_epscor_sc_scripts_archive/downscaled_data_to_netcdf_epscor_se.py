# convert the downscaled data archive 

def run( x ):
	''' simple wrapper to open and return a 2-D array from a geotiff '''
	import rasterio
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
# seasonal calculations
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

# def cf_attrs( scenario, model, contact='Michael Lindgren - malindgren@alaska.edu',  ):
# 	'''
# 	generate the cf_metadata convention attributes for the NC file
	
# 	CONVENTION SPEC HERE:
# 		http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html

# 	'''
# 	{'institution': 'Scenarios Network for Alaska + Arctic Planning' ,
# 	'institute_id':  'SNAP',
# 	'experiment_id':scenario,
# 	'source':model,
# 	'model_id':model,
# 	'forcing':,
# 	'parent_experiment_id':  ,
# 	'parent_experiment_rip':  ,
# 	'branch_time':  ,
# 	'contact':contact,
# 	'references': ,
# 	'initialization_method':  ,
# 	'physics_version':  ,
# 	'tracking_id':  ,
# 	'acknowledgements':  ,
# 	'cesm_casename':  ,
# 	'cesm_repotag':  ,
# 	'cesm_compset':  ,
# 	'resolution':  ,
# 	'forcing_note':  ,
# 	'processed_by':  ,
# 	'processing_code_information':  ,
# 	'product':  ,
# 	'experiment':  ,
# 	'frequency':  ,
# 	'creation_date':  ,
# 	'history':  ,
# 	'Conventions':'CF-1.6' ,
# 	'project_id':  ,
# 	'table_id':  ,
# 	'title':  ,
# 	'parent_experiment':  ,
# 	'modeling_realm':  ,
# 	'realization':  ,
# 	'cmor_version':  }

def generate_nc( model, variable, scenario, base_path, output_base_path, begin, end ):
	'''
	main function to output a netcdf file from a group of
	GeoTiff files of downscaled SNAP data.

	[MORE DOCS TO COME]
	'''
	# from pathos.multiprocessing import Pool
	from multiprocessing import Pool
	import numpy as np
	import pandas as pd
	import os, glob, rasterio, time, itertools
	import xarray as xr
	
	print( 'working on: {} {} {}'.format( variable, model, scenario ) )
	
	# set up pathing
	input_path = os.path.join( base_path, model, scenario, variable )
	output_path = os.path.join( output_base_path, model, scenario, variable )

	try: # try:except to overcome some multiprocessing collision issues
		if not os.path.exists( output_path ):
			os.makedirs( output_path )
	except:
		pass

	# list the data
	l = sort_files( glob.glob( os.path.join( input_path, '*.tif' ) ) )
	l = only_years( l, begin=begin, end=end )
	
	# open a pool and turn the list of arrays into an ndarray
	pool = Pool( ncpus )
	arr = np.array( pool.map( run, l ) )
	pool.close()
	pool.join()

	# mask it
	arr = np.ma.masked_where( arr <= np.min( arr ), arr )

	# [RECENT ADDITION] swap the axes so we are (lat, lon, time)
	arr = np.swapaxes( np.swapaxes(arr, 0, 2), 0, 1)

	# get the lons and lats for the NetCDF
	lons, lats = coordinates( l[0] )
	rst = rasterio.open( l[0] )

	# THIS IS A TEST AREA FOR PRODUCING THE *_bnds variables -- NOT IMPLEMENTED
	# # the res is standard in both directions.
	# res = 2000.0
	# half_res = 2000.0 / 2
	# lon_bnds = [ [i-half_res,i+half_res ] for i in lons.ravel() ]
	# # the lat_bnds variable appears to be the same as the above, but it is
	# # forced to the extent of the map so the lat_bnds at the top and bottom are 
	# # different resolution (half) of the remainder of the rectilinear grid cells.
	# # this needs to be taken into account in this calculation.
	# # MAYBE JUST HOLD IT TO THE EXTENT FOR THESE LATITUDES?
	# lat_bnds = [ [i-half_res,i+half_res ] for i in lats.ravel() ]
	# lat_mins, lat_max = rst.bounds 

	# get some time and date stuff
	t = time.time()

	# OGC WKT for EPSG:3338 which is the CF standard.
	crs_wkt = 'PROJCS["NAD83 / Alaska Albers",GEOGCS["NAD83",DATUM["North_American_Datum_1983",\
				SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],\
				AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],\
				AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",55],\
				PARAMETER["standard_parallel_2",65],PARAMETER["latitude_of_center",50],PARAMETER["longitude_of_center",-154],\
				PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],\
				AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3338"]]'

	# create the dataset in xarray
	ds = xr.Dataset( { variable:(['x','y','time'], arr) }, 
				coords={ 'lon': (['x', 'y'], lons),
						'lat': (['x', 'y'], lats),
						'time': pd.date_range( str(begin), str(end + 1), freq='M' ) },
				attrs={ 'units':'Celcius', 'time_interval':'monthly', 
						'variable':variable, 'model':model, 'scenario':scenario, 
						'crs_wkt':crs_wkt,
						'creation_date':time.ctime( t ), 'creation_date_UTC':t, 
						'created by':'Michael Lindgren - malindgren@alaska.edu',
						'nodata_value':'-3.39999995e+38',
						'cell_resolution':'2000 meters' } )

	# write it out to disk
	encoding = { variable: { '_FillValue':-3.39999995e+38, 'zlib':True } }
	output_filename = os.path.join( output_path, '_'.join([ variable, model, scenario, str( begin ), str( end ) ]) + '.nc' )
	ds.to_netcdf( output_filename, mode='w', encoding=encoding )
	ds.close() # close it
	return output_filename


if __name__ == '__main__':
	import os, glob
	import argparse

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="cmip5 model name (exact)" )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="cmip5 variable name (exact)" )
	parser.add_argument( "-s", "--scenario", action='store', dest='scenario', type=str, help="cmip5 scenario name (exact)" )
	args = parser.parse_args()

	# setup args
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5'
	output_base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_netcdf'
	units = 'C'
	time_interval = 'monthly'	
	ncpus = 32

	if args.scenario == 'historical':
		begin = 1900
		end = 2005
	else:
		begin = 2006
		end = 2100

	# main
	_ = generate_nc( args.model, args.variable, args.scenario, base_path, output_base_path, begin, end )
	