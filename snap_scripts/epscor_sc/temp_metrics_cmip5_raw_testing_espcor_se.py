# # # # # select from NetCDF Dataset across an AOI # # # # 

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
	of the month, which makes thing go 1, 11, ... which stinks for timeseries
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
	base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped'
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/testing_compare_raw'
	os.chdir( base_dir )
	model = 'MRI-CGCM3' #'IPSL-CM5A-LR' # 
	scenario = 'rcp45'
	variables = ['tasmax', 'tasmin', 'tas']
	shp_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/Kenai_StudyArea.shp'
	begin = 2030
	end = 2040

	# # shapefile work:
	# # warp to wgs84 greenwich, then to wgs84 pacific 
	shp = gpd.read_file( shp_fn )
	shp_ll = shp.to_crs( epsg=4326 )
	shp_pcll = shp_ll.to_crs( crs={ 'proj':'longlat','ellps':'WGS84','pm':-180,'datum':'WGS84' } )
	shapes = zip( shp_pcll.geometry, shp_pcll.OBJECTID )

	out = []
	for variable in variables:
		fn, = glob.glob( os.path.join( base_dir, model, scenario, variable, '*.nc' ) )
		ds = xr.open_dataset( fn ).load() - 273.15
		# mask it
		ds[ 'mask' ] = rasterize( shapes, ds.coords, longitude='lon', latitude='lat', fill=0, all_touched=True )
		# get the mean across space within the AOI for all timesteps
		mask_mean = ds[variable].sel( time=slice( str(begin), str(end) ) ).where( ds.mask == 1 ).mean( axis=(1,2) )	
		df = pd.DataFrame( mask_mean.to_pandas(), columns=[variable] )
		df.columns = [ variable ]
		out = out + [ df ]
		
	plot_df = pd.concat( out, axis=1 )
	plot_df = plot_df[['tasmax', 'tas', 'tasmin']]

	# now plot the dataframe
	if begin == end:
		title = 'EPSCoR SC AOI Temp Metrics raw {} {} {}'.format( model, scenario, begin )
	else:
		title = 'EPSCoR SC AOI Temp Metrics raw {} {} {} - {}'.format( model, scenario, begin, end )

	figsize = (13,9)
	colors = ['red', 'black', 'blue']

	ax = plot_df.plot( kind='line', title=title, figsize=figsize, color=colors )
	
	# now plot the dataframe
	if begin == end:
		plt.savefig( os.path.join( output_path,'mean_temps_epscor_sc_raw_{}_{}_{}.png'.format( model, scenario, begin ) ), dpi=600 )
	else:
		plt.savefig( os.path.join( output_path,'mean_temps_epscor_sc_raw_{}_{}_{}_{}.png'.format( model, scenario, begin, end ) ), dpi=600 )
	
	plt.close()


