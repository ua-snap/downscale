# # # # # compare tasmin, tas, tasmax in a timeseries of GeoTiff files # # # # 

def transform_from_latlon( lat, lon ):
	''' simple way to make an affine transform from lats and lons coords '''
	from affine import Affine
	lat = np.asarray( lat )
	lon = np.asarray(lon)
	trans = Affine.translation(lon[0], lat[0])
	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
	return trans * scale
def rasterize( shapes, coords, latitude='latitude', longitude='longitude', fill=None, **kwargs ):
	'''
	Rasterize a list of (geometry, fill_value) tuples onto the given
	xarray coordinates. This only works for 1d latitude and longitude
	arrays.
	'''
	from rasterio import features
	if fill == None:
		fill = np.nan
	transform = transform_from_latlon( coords[ latitude ], coords[ longitude ] )
	out_shape = ( len( coords[ latitude ] ), len( coords[ longitude ] ) )
	raster = features.rasterize(shapes, out_shape=out_shape,
								fill=fill, transform=transform,
								dtype=float, **kwargs)
	spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
	return xr.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))

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
	months = [ int(os.path.basename( fn ).split('.')[0].split( split_on )[elem_month]) for fn in files ]
	years = [ int(os.path.basename( fn ).split('.')[0].split( split_on )[elem_year]) for fn in files ]
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
	years = [ int(os.path.basename( fn ).split('.')[0].split( split_on )[elem_year]) for fn in files ]
	df = pd.DataFrame( { 'fn':files, 'year':years } )
	df_slice = df[ (df.year >= begin ) & (df.year <= end ) ]
	return df_slice.fn.tolist()
def masked_mean( fn ):
	''' get mean of the full domain since the data are already clipped 
	mostly used for processing lots of files in parallel.'''
	import numpy as np
	import rasterio

	with rasterio.open( fn ) as rst:
		mask = (rst.read_masks( 1 ) == 0)
		arr = np.ma.masked_array( rst.read( 1 ), mask=mask )
	return np.mean( arr )

if __name__ == '__main__':
	import os, glob
	import geopandas as gpd
	import numpy as np
	import xarray as xr
	import matplotlib
	matplotlib.use( 'agg' )
	from matplotlib import pyplot as plt
	from pathos.mp_map import mp_map
	import pandas as pd
	
	# args / set working dir
	base_dir = '/Users/malindgren/Documents/downscale_epscor/august_fix'
	os.chdir( base_dir )
	model = '5ModelAvg'
	scenario = 'rcp45'
	begin = 2010
	end = 2015

	variables = ['tasmax', 'tas', 'tasmin']
	out = {}
	for v in variables:
		path = os.path.join( base_dir,'EPSCOR_SC_DELIVERY_AUG2016','downscaled', model, scenario, v )
		# for testing with new downscaler
		if v == 'tas':
			path = os.path.join( base_dir,'downscaled_tas_pr_epscor_sc', model, scenario, v )
		
		files = glob.glob( os.path.join( path, '*.tif' ) )
		files = sort_files( only_years( files, begin=begin, end=end, split_on='_', elem_year=-1 ) )
		out[ v ] = mp_map( masked_mean, files, nproc=4 )

	plot_df = pd.DataFrame( out )
	plot_df.index = pd.date_range( start=str(begin), end=str(end+1), freq='M' )
	plot_df = plot_df[['tasmax', 'tas', 'tasmin']] # get em in the order for plotting

	# now plot the dataframe
	if begin == end:
		title = 'EPSCoR SC AOI Temp Metrics {} {} {}'.format( model, scenario, begin )
	else:
		title = 'EPSCoR SC AOI Temp Metrics {} {} {} - {}'.format( model, scenario, begin, end )

	figsize = (13,9)
	colors = ['red', 'black', 'blue' ]

	ax = plot_df.plot( kind='line', title=title, figsize=figsize, color=colors )
	# now plot the dataframe
	if begin == end:
		plt.savefig( 'mean_temps_epscor_sc_{}_{}_{}.png'.format( model, scenario, begin ), dpi=600 )
	else:
		plt.savefig( 'mean_temps_epscor_sc_{}_{}_{}_{}.png'.format( model, scenario, begin, end ), dpi=600 )
	
	plt.close()


# # # PRISM TEST VERSION DIFFERENCES # # # # # # #
# import rasterio
# import numpy as np
# import os, glob, itertools

# base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/raw_prism'
# variables = [ 'tmax', 'tmin' ]

# for variable in variables:
# 	ak_olds = sorted( glob.glob( os.path.join( base_path, 'prism_raw_older', 'ak', variable, '*.asc' ) ) )
# 	ak_news = sorted( glob.glob( os.path.join( base_path, 'prism_raw_2016', 'ak', variable, '*.asc' ) ) )

# 	olds = np.array([ rasterio.open( i ).read( 1 ) for i in ak_olds if '_14' not in i ])
# 	news = np.array([ rasterio.open( i ).read( 1 ) *.10 for i in ak_news if '_14' not in i ])

# 	out = olds - news
# 	out[ (olds == -9999.0) | (news == -9999.0) ] = 0

# 	uniques = np.unique( out )
# 	uniques[ uniques > 0.01 ]
