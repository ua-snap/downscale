# # # # # select from NetCDF Dataset across an AOI # # # # 

def transform_from_latlon( lat, lon ):
	''' simple way to make an affine transform from lats and lons coords '''
	from affine import Affine
	lat = np.asarray( lat )
	lon = np.asarray(lon)
	trans = Affine.translation(lon[0], lat[0])
	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
	return trans * scale
# def rasterize( shapes, coords, latitude='latitude', longitude='longitude', fill=None, **kwargs ):
# 	'''
# 	Rasterize a list of (geometry, fill_value) tuples onto the given
# 	xarray coordinates. This only works for 1d latitude and longitude
# 	arrays.
# 	'''
# 	from rasterio import features
# 	if fill == None:
# 		fill = np.nan
# 	transform = transform_from_latlon( coords[ latitude ], coords[ longitude ] )
# 	out_shape = ( len( coords[ latitude ] ), len( coords[ longitude ] ) )
# 	raster = features.rasterize(shapes, out_shape=out_shape,
# 								fill=fill, transform=transform,
# 								dtype=float, **kwargs)
# 	spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
# 	return xr.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))

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
	from rasterio.features import rasterize
	
	# args / set working dir
	# base_dir = '/Users/malindgren/Documents/downscale_epscor/august_fix'
	# os.chdir( base_dir )

	nc_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5_nwt_v2/cmip5_raw_ncrcat/GFDL-CM3/rcp45/tas/tas_Amon_GFDL-CM3_rcp45_200601-230012.nc'
	# nc_fn_hist = './tasmax_CCSM4_historical_r1i1p1_1850_2005.nc'
	shp_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/TEST_NWT/extent_rough_PCLL.shp'
	shp = gpd.read_file( shp_fn )
	begin = 2010
	end = 2019
	climbegin = 1961
	climend = 1990
	variable = 'tas'

	# read in data and convert to Celcius
	ds = xr.open_dataset( nc_fn ) #.load() - 273.15
	dat = ds[ variable ].data
	
	# make a gd geotransform from the input infos
	transform = transform_from_latlon( ds.lat.data.tolist(), ds.lon.data.tolist() )

	# make a raster mask with the shapefile using a single slice of the ds
	template_arr = ds[ variable ].isel(time=0).data

	# mask it 
	shapes = zip( shp.geometry, shp.id )
	mask = rasterize( shapes, out_shape=template_arr.shape, fill=0, transform=transform, dtype='float32' )
	mask_b = np.broadcast_to( mask, shape=dat.shape )

	# get the means through time
	dat_masked = np.ma.masked_where( mask_b == 1, dat, copy=True )
	monmeans = np.mean( dat_masked, axis=(1,2) ) - 273.15

	# make a df with the months, years, decades, and the monmeans
	years = range( 2006, 2300+1 )
	months = range( 1, 12+1 )

	df = pd.DataFrame([ {'month':m, 'year':y} for y in years for m in months ])
	df[ 'mean' ] = monmeans
	df[ 'decade' ] = df.year.apply( lambda x: int(str(x)[:3]+'0' ) )


	# warp to wgs84 greenwich, then to wgs84 pacific 
	# shp_ll = shp.to_crs( epsg=4326 )
	# shp_pcll = shp_ll.to_crs( crs={ 'proj':'longlat','ellps':'WGS84','pm':-180,'datum':'WGS84' } )
	shapes = zip( shp.geometry, shp.id )
	ds[ 'mask' ] = rasterize( shapes, ds.coords, longitude='lon', latitude='lat', fill=0 )

	# get the mean across space within the AOI for all timesteps
	orig_mean = ds.tas.where( ds.mask == 1 ).mean( axis=(1,2) )
	orig_mean = orig_mean - 273.15
	dat = orig_mean.data
	years = range(2006, 2300+1)
	# make a df and use it to group across space
	df = pd.DataFrame([ {'month':m,'year':y} for y in years for m in range(1,12+1) ])
	df['dat'] = dat
	df.decade = df.year.apply(lambda x: int(str(x)[:3]+'0'))

	df.groupby(df.decade).apply(np.mean, axis=0)



	# # # # # CALCULATE ANOMALIES FROM HISTORICAL of same model/scenario/variable
	ds_hist[ 'mask' ] = rasterize( shapes, ds_hist.coords, longitude='lon', latitude='lat', fill=0 )
	
	# climatology = ds_hist.tasmax.sel( time=slice( str(climbegin), str(climend) ) ).groupby( 'time.month' ).mean( axis=0 )
	climatology = ds_hist.tasmax.sel( time=slice( str(climbegin), str(climend) ) )
	climatology = climatology.groupby( 'time.month' ).apply( lambda x: np.mean( x, axis=0 ) )
	anomalies = ds.tasmax.groupby( 'time.month' ) - climatology
	anomalies[ 'mask' ] = rasterize( shapes, anomalies.coords, longitude='lon', latitude='lat', fill=0 )
	anom_mean = anomalies.sel( time=slice( str(begin), str(end) ) ).where( anomalies.mask == 1 ).mean( axis=(1,2) )

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# now lets do the same thing with the rasterio-downscaled data
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	path = '/Users/malindgren/Documents/downscale_epscor/august_fix/EPSCOR_SE_DELIVERY_AUG2016/downscaled/5ModelAvg/rcp45/tasmax'
	# path = '/Users/malindgren/Documents/downscale_epscor/august_fix/EPSCOR_SC_DELIVERY_AUG2016/derived/grids/annual_seasonals/5ModelAvg/rcp45/tasmax'
	files = glob.glob( os.path.join( path, '*.tif' ) )
	files = sort_files( only_years( files, begin=begin, end=end, split_on='_', elem_year=-1 ) )
	
	# get the means across space for each timestep
	down_mean = mp_map( masked_mean, files, nproc=4 )
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# # # # # NOW LETS SEE THIS WITH THE LINEAR INTERPOLATION USED... as a TEST. # # # # # # # #
	# path = '/Users/malindgren/Documents/downscale_epscor/august_fix/CCSM4_clip/tasmax'
	# files = glob.glob( os.path.join( path, '*.tif' ) )
	# files = sort_files( only_years( files, begin=begin, end=end, split_on='_', elem_year=-1 ) )
	
	# # get the means across space for each timestep
	# down_mean_linear = mp_map( masked_mean, files, nproc=4 )
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# now lets do the same thing with the PRISM data
	# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
	# # clipped data with the below with the epscor crop/clip code:
	# base_path = '/Users/malindgren/Documents/downscale_epscor/tasmax'
	# output_path = '/Users/malindgren/Documents/downscale_epscor/tasmax_prism_epscor'
	# ncpus = 4
	# subdomain_fn = '/Users/malindgren/Documents/downscale_epscor/august_fix/Kenai_StudyArea.shp'
	path = '/Users/malindgren/Documents/downscale_epscor/tasmax_prism_epscor'
	files = glob.glob( os.path.join( path, '*.tif' ) )
	files = sorted( [ i for i in files if '_14' not in i ] )
	prism_mean = mp_map( masked_mean, files, nproc=4 )

	# NOW LETS PLOT BOTH OF THESE GUYS
	# make a pandas dataframe with the xarray dataset
	plot_df = pd.DataFrame( orig_mean.to_pandas() )
	plot_df.columns = ['CCSM4_prepped']
	plot_df['5ModelAvg_down'] = down_mean
	# plot_df['down_linear'] = down_mean_linear
	# plot_df['prism'] = prism_mean * (len( down_mean ) / 12)
	plot_df['CCSM4_anom'] = anom_mean

	# now plot the dataframe
	plot_df.plot( kind='line' )
	plt.savefig( 'mean_southcentral_testcase_fixed_5model.png' )
	plt.close()



# /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prepped_cmip5/CCSM4/rcp45/tasmax/tasmax_CCSM4_rcp45_r1i1p1_2006_2100.nc
# /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/Kenai_StudyArea.*

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # NOTES AND STUFF # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # QUESTIONS FOR MATT:

# 1. How we calculate the metrics
# 	- currently I am calculating the mean through time, then a metric across the resulting
# 		2D array.

# 2. What to do with the NoData areas in the Northern portion of the AOI, which is an artifact
# 	of proportional downscaling.
# 	- maybe just bring PRISM through those locations?
# 	- seems like doing the above will cause for some linear anomalies in the output raster from the 
# 		low resolution 'holes' in the landscape.

# 	** - check the INPUT data BEFORE the division.  This needs to be checked. 0.5 fill-in values. - **

# 3. What do you think we should do to examine these extreme outputs?  John Walsh and some other people
# 	who I do not know, seem to think that the data is 'not right', of course without ever getting 
# 	their hands dirty they say that the data doesnt seem to mesh well with the mean data.  
# 	- I have already dug into the code to make sure it is grabbing the files correctly, and it 
# 		appears to be working correctly.  There is a step where I have to forcibly make the files into
# 		a single array with the _right_ years and such.  This standardization is important for working 
# 		with the data in the varied NetCDF temporal domains that each of the models slice the 
# 		data into.
# 	- 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# FROM MATT:
# ----------
# -> GCM anomalies -- add to PRISM -- get the domain in the GCMs.
# -> AKCAN full domain.  Color gradient of the quantity of unreasonable values.
# -> 
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


