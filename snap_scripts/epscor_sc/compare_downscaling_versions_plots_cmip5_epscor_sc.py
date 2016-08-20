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
def masked_mean( fn, bounds=None ):
	''' get mean of the full domain since the data are already clipped 
	mostly used for processing lots of files in parallel.'''
	import numpy as np
	import rasterio
		
	with rasterio.open( fn ) as rst:
		if bounds:
			window = rst.window( *bounds )
		else:
			window = rst.window( *rst.bounds )
		mask = (rst.read_masks( 1 ) == 0)
		arr = np.ma.masked_array( rst.read( 1, window=window ), mask=mask )
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
	import geopandas as gpd
	
	# args / set working dir
	base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data'
	os.chdir( base_dir )
	scenario = 'historical'
	shp_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/Kenai_StudyArea.shp'
	shp = gpd.read_file( shp_fn )
	bounds = shp.bounds

	models = ['5ModelAvg','CRU_TS323','GFDL-CM3','GISS-E2-R','IPSL-CM5A-LR','MRI-CGCM3','NCAR-CCSM4']
	variables_list = [['tasmax', 'tas', 'tasmin'],['pr']]

	for variables in variables_list:
		for m in models:
			if m == 'CRU_TS323':
				old_dir = '/Data/Base_Data/Climate/AK_CAN_2km/historical/CRU/CRU_TS32'
				begin = 2000
				end = 2010
			else:
				if scenario == 'historical':
					old_dir = '/Data/Base_Data/Climate/AK_CAN_2km/historical/AR5_CMIP5_models'
					begin = 2000
					end = 2005
				else:
					old_dir = '/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models'
					begin = 2010
					end = 2015

			figsize = (16,9)
			out = {}
			for v in variables:
				path = os.path.join( base_dir,'downscaled', m, scenario, v )
				files = glob.glob( os.path.join( path, '*.tif' ) )
				files = sort_files( only_years( files, begin=begin, end=end, split_on='_', elem_year=-1 ) )
				out[ v ] = mp_map( masked_mean, files, nproc=4 )
				if v == 'tas' or v == 'pr':
					if m == 'CRU_TS323':
						path = os.path.join( old_dir, v )
					else:	
						path = os.path.join( old_dir, scenario, m, v )
					files = glob.glob( os.path.join( path, '*.tif' ) )
					files = sort_files( only_years( files, begin=begin, end=end, split_on='_', elem_year=-1 ) )
					out[ v+'_old' ] = mp_map( masked_mean, files, nproc=4 )

			plot_df = pd.DataFrame( out )
			plot_df.index = pd.date_range( start=str(begin), end=str(end+1), freq='M' )
			
			# sort the columns for output plotting cleanliness:
			if 'tas' in variables:
				col_list = ['tasmax', 'tas_old', 'tas', 'tasmin']
			elif 'pr' in variables:
				col_list = ['pr', 'pr_old']
			
			plot_df = plot_df[ col_list ] # get em in the order for plotting

			# now plot the dataframe
			if begin == end:
				title = 'EPSCoR SC AOI Temp Metrics {} {} {}'.format( m, scenario, begin )
			else:
				title = 'EPSCoR SC AOI Temp Metrics {} {} {} - {}'.format( m, scenario, begin, end )

			if v == 'tas':
				colors = ['darkred', 'black', 'blue', 'lightred' ]
			else:
				colors = [ 'blue', 'black' ]

			ax = plot_df.plot( kind='line', title=title, figsize=figsize, color=colors )

			output_dir = os.path.join( base_dir, 'compare_downscaling_vars' )
			if not os.path.exists( output_dir ):
				os.makedirs( output_dir )

			# now plot the dataframe
			out_metric_fn = 'temps'
			if 'pr' in variables:
				out_metric_fn = 'prec'

			if begin == end:
				output_filename = os.path.join( output_dir,'mean_{}_epscor_sc_{}_{}_{}.png'.format( out_metric_fn, m, scenario, begin ) )
			else:
				output_filename = os.path.join( output_dir,'mean_{}_epscor_sc_{}_{}_{}_{}.png'.format( out_metric_fn, m, scenario, begin, end ) )
			plt.savefig( output_filename, dpi=700 )
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
