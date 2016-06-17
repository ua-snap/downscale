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


if __name__ == '__main__':
	import glob, os, rasterio, itertools
	import numpy as np
	from pathos import multiprocessing
	import pandas as pd
	import xarray as xr

	# some setup args
	base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5'
	variables = [ 'tasmin', 'tasmax' ]
	scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ]

	for variable, model, scenario in itertools.product( variables, models, scenarios ):
		seasons = {'DJF':[12, 1, 2],'MAM':[3, 4, 5],'JJA':[6, 7, 8],'SON':[9, 10, 11]}
		# list the files
		files = glob.glob( os.path.join( base_dir, model, scenario, variable, '*.tif' ) )
		months = [ int(f.split('.')[0].split('_')[-2]) for f in files ]
		years = [ int(f.split('.')[0].split('_')[-1]) for f in files ]
		df = pd.DataFrame( {'fn':files, 'month':months, 'year':years} )
		df = df.sort_values(['year', 'month'])
		
		# some begin end stuff
		begin = '-'.join( df.iloc[0].ix[['month','year']].astype( str ) )
		end = '-'.join( ['1', (df.iloc[-1].ix['year']+1).astype( str )] )

		# get the lons and lats for the NetCDF
		lons, lats = coordinates( files[0] )
		lons, lats = [ np.unique( x ) for x in [lons, lats] ]

		# make the huge data array
		pool = multiprocessing.Pool( 32 )
		arr = np.array( pool.map( lambda x: rasterio.open( x ).read( 1 ), files ) )
		pool.close()
		pool.join()
		break

		# create the dataset in xarray
		ds = xr.Dataset( { variable:(['time','x','y'], arr) }, 
					coords={ 'lon': (['x', 'y'], lons),
							'lat': (['x', 'y'], lats),
							'time': pd.date_range( begin, end, freq='M' ) },
					attrs={ 'units':'Celcius', 'time_interval':'monthly', 
							'variable':variable, 'model':model, 'scenario':scenario } )

		# now create a seasonal average...
		ds_season = ds.groupby( 'time.season' ).mean( axis=0 )
		