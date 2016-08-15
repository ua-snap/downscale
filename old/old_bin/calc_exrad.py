# Python re-write of the R calcRa formula from Steph McAfee
def coordinates( fn, latlong=True ):
	'''
	take a raster file as input and return the centroid coords for each 
	of the grid cells as a pair of numpy 2d arrays (longitude, latitude)
	'''
	import rasterio
	import numpy as np
	from affine import Affine
	from pyproj import Proj, transform

	# Read raster
	with rasterio.open(fn) as r:
		T0 = r.affine  # upper-left pixel corner affine transform
		p1 = Proj(r.crs)
		A = r.read( 1 )  # pixel values

	# All rows and columns
	cols, rows = np.meshgrid( np.arange( A.shape[1] ), np.arange( A.shape[0] ) )

	# Get affine transform for pixel centres
	T1 = T0 * Affine.translation( 0.5, 0.5 )
	
	# Function to convert pixel row/column index (from 0) to easting/northing at centre
	rc2en = lambda r, c: (c, r) * T1

	# All eastings and northings (there is probably a faster way to do this)
	eastings, northings = np.vectorize(rc2en, otypes=[np.float, np.float])(rows, cols)

	if latlong != True:
		return eastings, northings
	elif latlong -- True:
		# Project all longitudes, latitudes to latlong
		longs, lats = transform(p1, p1.to_latlong(), eastings, northings)
		return longs, lats

def calcRa( lats, jd ):
	'''
	################################################################################
	This function calculates extraterrestrial solar radiation (radiation at the top
	of the earth's atmosphere, aka insolation) from day of the year and latitude.
	It is calculated separately for each day of the year.
	NB: Latitude should be in radians (not degrees).
	NB: Non-real values are produced for w (sunset hour angle) at times of the year
	and at latitudes where the sun never rises or sets. However, it makes sense to
	use just the real portion of the value, as it sets w to 0 in the winter and pi
	(~3.14) in the summer. 
	Source: Allen et al. (1998)

	>>> Translated to Python from Stephanie McAfee's R script by Michael Lindgren <<<
	#################################################################################

	Arguments:
		lats = {numpy.ndarray} latitudes of raster pixel centroids
		jd = {int} julian day (1 - 365)
	'''
	# convert lats to radians:
	lats = (lats*np.pi)/180

	# Calculate the earth-sun distance, which is a function solely of Julian day.  It is a single value for each day
	D = 1 + ( 0.033 * np.cos( 2 * np.pi * jd / 365 ) )

	# Calculate declination, a function of Julian day.  It is a single value for each day.
	DC = 0.409 * np.sin( ( ( 2 * np.pi / 365) * jd ) - 1.39 )

	# Calculate the sunset hour angle, a function of latitude and declination. Note that at v. high latitudes, this function can produce non-real values.  
	w = np.arccos( ( -1 * np.tan( DC ) * np.tan( lats ) ).astype( 'complex64' ) ).real
	
	# Calculate Ra using the above variables
	# S is the solar constant and is =0.082 MJ/m2min
	ra = (24 * 60 / np.pi) * D * 0.082 *( w * np.sin(lats) * np.sin(DC) + np.cos(lats) * np.cos(DC) * np.sin(w) )
	return ra

def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

if __name__ == '__main__':
	from collections import OrderedDict
	import rasterio, os
	import numpy as np
	import pandas as pd
	import geopandas as gpd
	from shapely.geometry import Point

	output_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/girr_radiation'
	template_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'

	# some setup of the grouping of the Julian Days to Months
	month_days_dict = OrderedDict( [('01',31),('02',28),('03',31),('04',30),('05',31),('06',30),('07',31),('08',31),('09',30),('10',31),('11',30),('12',31)] )
	julian_month_splitter = np.hstack(np.array([ np.repeat( count+1, i ).tolist() for count, i in enumerate( month_days_dict.values() ) ]))
	julian_days = pd.Series( range( 1, 366 ) )
	julian_months_grouped = julian_days.groupby( julian_month_splitter )

	template_rst = rasterio.open( template_fn )
	meta = template_rst.meta
	meta.update( compress='lzw' )
	lons, lats = coordinates( template_fn, False ) # test change
	lats = np.ma.masked_where( template_rst.read( 1 ) < -3.39999995e+34, lats )

	for month, days in julian_months_grouped.indices.iteritems():
		month = str( month )
		if len( month ) == 1:
			month = '0' + month
		month_mean = np.dstack([ calcRa( lats, day ) for day in days+1 ]).mean( axis=2 )
		
		# [TEST]: it may be necessary to take this array and its coords in latlong and reproject it to 3338
		lons = np.ma.masked_where( template_rst.read( 1 ) < -3.39999995e+34, lats )

		pts = [ Point( lalo ) for lalo in zip(lons.ravel().tolist(), lats.ravel().tolist()) ]
		mm = month_mean.ravel().tolist()
		df = pd.DataFrame( {'Ra':mm, 'geometry':pts} ) # remove masking
		gdf = gpd.GeoDataFrame( df )
		break

		# then rasterize this to the extent of the template_rst

		# [END TEST]

		output_filename = os.path.join( output_path, 'ra_mean_allen1998_'+month+'netest_akcan.tif' )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( month_mean.astype( template_rst.dtypes[ 0 ] ), 1 )


# # what does R say? <<- I think it is saying that we did it slightly wrong, but the above code can be used to generate proper outputs
# require( raster )

# template_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
# r <- as.matrix( coordinates( raster( template_fn ) ) )
# jd = 100

