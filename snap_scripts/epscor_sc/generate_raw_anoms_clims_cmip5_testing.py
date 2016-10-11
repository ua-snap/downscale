# WORK WITH MATT:
def transform_from_latlon( lat, lon ):
	''' simple way to make an affine transform from lats and lons coords '''
	from affine import Affine
	lat = np.asarray( lat )
	lon = np.asarray( lon )
	trans = Affine.translation(lon[0], lat[0])
	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
	return trans * scale

if __name__ == '__main__':
	import rasterio, os
	import xarray as xr
	import numpy as np
	import pandas as pd

	# some id variables
	variable = 'tas'
	model = 'GFDL-CM3'
	scenario = 'rcp60'
	begin = '1961'
	end = '1990'
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/matt_mike_test'

	# read em
	historical = xr.open_dataset( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/{}/{}/{}/{}_{}_{}_r1i1p1_1860_2005.nc'.format( model, 'historical', variable, variable, model, 'historical' ) )

	future = xr.open_dataset( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/{}/{}/{}/{}_{}_{}_r1i1p1_2006_2100.nc'.format( model, scenario, variable, variable, model, scenario ) )

	# concat em
	ds = xr.concat( [historical, future], dim='time' )

	# make climatology
	climatology = ds[variable].sel( time=slice( begin, end ) ).groupby( 'time.month' ).mean( 'time' )

	# make anomalies
	anomalies = ds[variable].groupby( 'time.month' ) - climatology

	# slice anomalies back to just the future
	anomalies = anomalies.sel( time=future.time )

	# make some metadata
	count, height, width = anomalies.shape
	affine = transform_from_latlon( ds.lat, ds.lon )
	meta = { 'crs':{'init':'epsg:4326'},
			 'count':count,
			 'height':height,
			 'width':width,
			 'driver':'GTiff',
			 'dtype':'float64',
			 'affine':affine }

	output_filename = os.path.join( output_path, '{}_{}_{}_multiband_anomalies.tif'.format( variable, model, scenario ) )

	# write anom to disk
	with rasterio.open( output_filename, 'w', **meta ) as rst:
		rst.write( anomalies.values )

	# write climatologies to disk
	# make some metadata
	count, height, width = climatology.shape
	affine = transform_from_latlon( ds.lat, ds.lon )
	meta = { 'crs':{'init':'epsg:4326'},
			 'count':count,
			 'height':height,
			 'width':width,
			 'driver':'GTiff',
			 'dtype':'float64',
			 'affine':affine }

	output_filename = os.path.join( output_path, '{}_{}_{}_multiband_climatology.tif'.format( variable, model, scenario ) )

	# write anom to disk
	with rasterio.open( output_filename, 'w', **meta ) as rst:
		rst.write( climatology.values )
