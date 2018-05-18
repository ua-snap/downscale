# # # # # # # # # # # # # # # # # # # # # # # #
# PORT S.McAffee's Ra SCRIPT TO Python
# # # # # # # # # # # # # # # # # # # # # # # #

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
			T0 = r.transform  # upper-left pixel corner affine transform
			p1 = Proj( r.crs )
			A = r.read( 1 )  # pixel values

	elif (meta is not None) & (numpy_array is not None):
		A = numpy_array
		if input_crs != None:
			p1 = Proj( input_crs )
			T0 = meta[ 'transform' ]
		else:
			p1 = None
			T0 = meta[ 'transform' ]
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

def calc_ra( day, lat ):
	'''
	calculate Ra (a direct port to Python from S.McAfee R script) based on Allen et.al 1998

	ARGUMENTS:
	----------
	day = [int] Ordinal (1-365*) Day of the year to compute Ra
	lat = [np.ndarray] 2-D Numpy array with Latitude values converted to radians

	RETURNS:
	--------
	numpy.ndarray containing Ra values over the AOI of `lat`

	'''
	import numpy as np
	
	#Calculate the earth-sun distance, which is a function solely of Julian day.  It is a single value for each day
	d = 1+(0.033*np.cos( (2*np.pi*day/365) ) )

	#Calculate declination, a function of Julian day.  It is a single value for each day.
	dc = 0.409*np.sin(((2*np.pi/365)*day)-1.39)
	w = np.nan_to_num(np.real( np.arccos( ( -1*np.tan( dc )*np.tan( lat ) ) ).astype(np.complex_)))
	return (24*60/np.pi) * d * 0.082 * (w*np.sin(lat)*np.sin(dc)+np.cos(lat)*np.cos(dc)*np.sin(w))


if __name__ == '__main__':
	import rasterio, datetime, os
	import numpy as np
	import pandas as pd
	import geopandas as gpd
	from functools import partial
	from pathos.mp_map import mp_map
	from shapely.geometry import Point
	from pyproj import Proj, transform

	fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/insolation_L48/prism_raw_template/PRISM_tmean_30yr_normal_800mM2_01_bil.tif'
	output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/insolation_L48/climatologies'
	lons, lats = coordinates( fn )

	if not os.path.exists( output_path ):
		os.makedirs( output_path )

	rst = rasterio.open( fn )

	# mask those lats so we dont compute where we dont need to:
	data_ind = np.where( rst.read_masks( 1 ) != 0 )
	pts = zip( lons[ data_ind ].ravel().tolist(), lats[ data_ind ].ravel().tolist() )

	# radians from pts
	p1 = Proj( init='epsg:4269' )
	p2 = Proj( init='epsg:4326' )
	transform_p = partial( transform, p1=p1, p2=p2 )
	pts_radians = [ transform_p( x=lon, y=lat, radians=True ) for lon,lat in pts ]
	lat_rad = pd.DataFrame( pts_radians, columns=['lon','lat']).lat

	# calc ordinal days to compute
	ordinal_days = range( 1, 365+1, 1 )
	# make a monthly grouper of ordinal days
	ordinal_to_months = [ str(datetime.date.fromordinal( i ).month) for i in ordinal_days ]
	# convert those months to strings
	ordinal_to_months = [ ('0'+month if len( month ) < 2 else month) for month in ordinal_to_months  ]

	# calc girr
	f = partial( calc_ra, lat=lat_rad )
	Ra = mp_map( f, ordinal_days, nproc=32 )
	Ra_monthlies = pd.Series( Ra ).groupby( ordinal_to_months ).apply( lambda x: np.array(x.tolist()).mean( axis=0 ) )

	# iteratively put them back in the indexed locations we took them from
	meta = rst.meta
	meta.update( compress='lzw', count=1, dtype='float32' )
	for month in Ra_monthlies.index:
		arr = rst.read( 1 )
		arr[ data_ind ] = Ra_monthlies.loc[ month ].tolist()
		output_filename = os.path.join( output_path, 'girr_w-m2_{}.tif'.format(str( month ) ) )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( arr.astype( np.float32 ), 1 )
