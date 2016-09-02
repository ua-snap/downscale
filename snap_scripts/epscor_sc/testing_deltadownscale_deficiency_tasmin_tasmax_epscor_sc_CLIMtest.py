# SLICE YEARS AND CROP TO THE BASE EXTENT OF AKCAN IN WGS84 PCLL
def transform_from_latlon( lat, lon ):
	''' simple way to make an affine transform from lats and lons coords '''
	from affine import Affine
	lat = np.asarray( lat )
	lon = np.asarray( lon )
	trans = Affine.translation(lon[0], lat[0])
	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
	return trans * scale

if __name__ == '__main__':
	import matplotlib
	matplotlib.use( 'agg' )
	from matplotlib import pyplot as plt
	import os, rasterio
	from rasterio import crs
	import xarray as xr
	import numpy as np
	import pandas as pd
	import geopandas as gpd

	filelist = [ '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/IPSL-CM5A-LR/rcp85/tasmin/tasmin_IPSL-CM5A-LR_rcp85_r1i1p1_2006_2100.nc',
				'/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/IPSL-CM5A-LR/rcp85/tasmax/tasmax_IPSL-CM5A-LR_rcp85_r1i1p1_2006_2100.nc', 
				'/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/IPSL-CM5A-LR/rcp85/tas/tas_IPSL-CM5A-LR_rcp85_r1i1p1_2006_2100.nc' ]

	historicals = [ '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/IPSL-CM5A-LR/historical/tasmin/tasmin_IPSL-CM5A-LR_historical_r1i1p1_1860_2005.nc',
				'/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/IPSL-CM5A-LR/historical/tasmax/tasmax_IPSL-CM5A-LR_historical_r1i1p1_1860_2005.nc', 
				'/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/IPSL-CM5A-LR/historical/tas/tas_IPSL-CM5A-LR_historical_r1i1p1_1860_2005.nc' ]

	variables = ['tasmin', 'tasmax', 'tas']
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_test_ipsl/IPSL-CM5A-LR_raw/rcp85'

	begin = '2096'
	end = '2096'
	model = 'IPSL-CM5A-LR'
	scenario = 'rcp85'
	figsize = (16,9)

	# pacific-centered reference system 4326
	pacific = "+proj=longlat +ellps=WGS84 +pm=-180 +datum=WGS84 +no_defs"
	shp_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/SCTC_studyarea/Kenai_StudyArea.shp'
	shp = gpd.read_file( shp_fn )
	shp_pacific = shp.to_crs( pacific ) 
	bounds = shp_pacific.bounds

	raw_dict = {}
	anom_dict = {}
	for variable, fn, fn_hist in zip( variables, filelist, historicals ):
		print( 'running: {}'.format( variable ) )
		ds = xr.open_dataset( fn )
		ds = ds[ variable ].sel( time=slice( begin, end ) )
		ds = ds - 273.15
		time_suffix = [ t.strftime('%m_%Y') for t in ds.time.to_pandas() ]

		# get the historicals for climatology
		ds_hist = xr.open_dataset( fn_hist )
		ds_hist = ds_hist[ variable ].sel( time=slice( '01-1961', '12-1990' ) )
		ds_hist = ds_hist - 273.15
		climatology = ds_hist.groupby( 'time.month' ).mean( axis=0 )

		# create anomalies
		anomalies = ds.groupby( 'time.month' ) - climatology

		count, height, width = climatology.shape
		affine = transform_from_latlon( ds.lat, ds.lon )

		# can we window the data here?
		# this will involve generating a window object using rasterio and extracting using slicing.
		meta = { 	'crs':crs.from_string( pacific ),
					'affine':affine,
					'dtype':'float64',
					'driver':'GTiff',
					'height': height,
					'width': width,
					'count':count	}
		
		# write out the climatlogies

		# # # # CLIMATOLOGY OUTPUT
		#out_fn = ''
		out_fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_test_ipsl/IPSL-CM5A-LR/'+ variable +'_IPSL-CM5A-LR_rcp85_clim.tif'
		with rasterio.open( out_fn, 'w', **meta ) as rst:
			window = rst.window( *bounds.as_matrix().ravel().tolist() )
			rst.write( climatology.data )

			# read it back from window
			new_affine = rst.window_transform( window )
			raw_arr = rst.read( window=window )

			# same for anom arr
			rst.write( anomalies )
			clim_arr = rst.read( window=window )


		out_fn = ''
		with rasterio.open( out_fn, 'w', **meta ) as rst:
			window = rst.window( *bounds.as_matrix().ravel().tolist() )
			rst.write( ds.data )

			# read it back from window
			new_affine = rst.window_transform( window )
			raw_arr = rst.read( window=window )

			# same for anom arr
			rst.write( anomalies )
			anom_arr = rst.read( window=window )

		rst = None

		count, height, width = raw_arr.shape
		meta.update( height=height, width=width, affine=new_affine, driver='GTiff' )

		# write out raw arr
		out_fn = os.path.join( output_path, '{}_IPSL_CM5A_raw_data_{}_epscor_sc_{}_{}.tif'.format( variable,'rcp85',begin, end ))
		dirname = os.path.dirname( out_fn )
		if not os.path.exists( dirname ):
			os.makedirs( dirname )

		with rasterio.open( out_fn, 'w', **meta ) as out:
			out.write( raw_arr )

		# write out anom_arr
		out_fn = os.path.join( output_path, '{}_IPSL_CM5A_raw_anom_data_{}_epscor_sc_{}_{}.tif'.format( variable,'rcp85',begin, end ))
		dirname = os.path.dirname( out_fn )
		if not os.path.exists( dirname ):
			os.makedirs( dirname )

		with rasterio.open( out_fn, 'w', **meta ) as out:
			out.write( anom_arr )

		raw_dict[ variable ] = raw_arr.mean( axis=(1,2) )
		anom_dict[ variable ] = anom_arr.mean( axis=(1,2) )

	# plot the above data in groups showing the raw values and the 
	df_raw = pd.DataFrame( raw_dict )
	df_anom = pd.DataFrame( anom_dict )
	
	col_list = ['tasmax', 'tas', 'tasmin']
	df_raw = df_raw[ col_list ]
	df_raw.index = pd.date_range('2096', '2097', freq='M')
	df_anom = df_anom[ col_list ]
	df_anom.index = pd.date_range('2096', '2097', freq='M')

	# RAW FIRST
	# now plot the dataframe
	if begin == end:
		title = 'EPSCoR SC AOI Temp Mean {} {} {}'.format( model, scenario, begin )
	else:
		title = 'EPSCoR SC AOI Temp Mean {} {} {} - {}'.format( model, scenario, begin, end )

	if 'tas' in variables:
		colors = ['red', 'black', 'blue' ]
	else:
		colors = [ 'blue', 'black' ]

	ax = df_raw.plot( kind='line', title=title, figsize=figsize, color=colors )

	output_dir = os.path.join( output_path, 'test_min_max_issue' )
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )

	out_metric_fn = 'temps'
	if 'pr' in variables:
		out_metric_fn = 'prec'

	if begin == end:
		output_filename = os.path.join( output_dir,'mean_{}_epscor_sc_{}_{}_{}.png'.format( out_metric_fn, model, scenario, begin ) )
	else:
		output_filename = os.path.join( output_dir,'mean_{}_epscor_sc_{}_{}_{}_{}.png'.format( out_metric_fn, model, scenario, begin, end ) )
	plt.savefig( output_filename, dpi=700 )
	plt.close()

	# ANOM NEXT
	# now plot the dataframe
	if begin == end:
		title = 'EPSCoR SC AOI Temp Mean Anomalies {} {} {}'.format( model, scenario, begin )
	else:
		title = 'EPSCoR SC AOI Temp Mean Anomalies {} {} {} - {}'.format( model, scenario, begin, end )

	if 'tas' in variables:
		colors = ['red', 'black', 'blue' ]
	else:
		colors = [ 'blue', 'black' ]

	ax = df_anom.plot( kind='line', title=title, figsize=figsize, color=colors )

	output_dir = os.path.join( output_path, 'test_min_max_issue' )
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )

	# now plot the dataframe
	out_metric_fn = 'temps'
	if 'pr' in variables:
		out_metric_fn = 'prec'

	if begin == end:
		output_filename = os.path.join( output_dir,'mean_{}_epscor_sc_anom_{}_{}_{}.png'.format( out_metric_fn, model, scenario, begin ) )
	else:
		output_filename = os.path.join( output_dir,'mean_{}_epscor_sc_anom_{}_{}_{}_{}.png'.format( out_metric_fn, model, scenario, begin, end ) )
	plt.savefig( output_filename, dpi=700 )
	plt.close()


