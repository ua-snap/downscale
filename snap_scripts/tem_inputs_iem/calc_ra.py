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
	import rasterio, datetime
	import numpy as np
	import pandas as pd
	import geopandas as gpd
	from functools import partial
	from pathos.mp_map import mp_map
	from shapely.geometry import Point

	fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/downscaled/CCSM4/rcp26/tas/tas_mean_C_ar5_NCAR-CCSM4_rcp26_01_2006.tif'
	lons, lats = coordinates( fn )

	rst = rasterio.open( fn )

	# mask those lats so we dont compute where we dont need to:
	# lons[ rst.read_masks( 1 ) == 0  ] = -9999
	lats[ rst.read_masks( 1 ) == 0  ] = -9999
	data_ind = np.where( lats != -9999 )

	# not yet working below:
	pts = pd.DataFrame( { 'lat':lats.ravel(), 'lon':lons.ravel() } )
	pts[ 'geometry' ] = gpd.GeoSeries( pts.apply( lambda x: Point(x.lon, x.lat), axis=1 ) )
	pts = gpd.GeoDataFrame( pts, crs={'init':'epsg:3338'}, geometry='geometry' )
	pts_ll = pts.to_crs( epsg=4326 )
	lat_rad = np.radians( pts_ll['lat'].reshape( *lons.shape ) )

	# forget the above for testing, lets use Stephs radians
	# latr = rasterio.open('/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/radiance/radians.txt')
	# latr = latr.read( 1 )

	# calc ordinal days to compute
	ordinal_days = range( 1, 365+1, 1 )
	ordinal_to_months = [ datetime.date.fromordinal( i ).month for i in ordinal_days ]

	f = partial( calc_ra, lat=lat_rad )
	Ra = mp_map( f, ordinal_days, nproc=32 )
	Ra_monthlies = pd.Series( Ra ).groupby( ordinal_to_months ).apply( lambda x: np.array(x.tolist()).mean( axis=0 ) )




# JUNK FOR NOW
# SOME SETUP
# # this little bit just makes some sample grids to use in calculations
# fn = '/Users/malindgren/Documents/downscale_epscor/sept_fix/CalcRa/test_calcRa_4326_small.tif'
# rst = rasterio.open( fn  )
# meta = rst.meta
# meta.update( compress='lzw', count=1, nodata=None )
# meta.pop( 'transform' )
# new_fn = fn.replace( '_small', '' )
# with rasterio.open( new_fn, 'w', **meta ) as out:
# 	out.write( rst.read(1), 1 )
# # # # # 