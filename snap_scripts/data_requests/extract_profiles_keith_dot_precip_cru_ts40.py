# extraction for Keith Cunningham -- Glitter Gulch DOT
def make_gdf():
	df = pd.DataFrame({'name':['Long Lake','Glitter Gulch'],
		'lat':[61.8,63.76],'lon':[-148.2,-148.9]})

	df['geometry'] = df.apply(lambda x: Point(x.lon, x.lat), axis=1)
	shp = gpd.GeoDataFrame( df, crs={'init':'epsg:4326'}, geometry='geometry')
	return shp.to_crs( epsg=3338 )

def list_data( base_dir ):
	files = glob.glob( os.path.join( cur_path, '*.tif' ) )
	df = pd.DataFrame([ os.path.basename(fn).split('.')[0].split('_')[-2:] for fn in files ], columns=['month','year'])
	df['fn'] = files
	return df.sort_values(['year','month']).reset_index(drop=True)

def rasterize( shapes, coords, latitude='lat', longitude='lon', fill=None, **kwargs ):
	'''
	Rasterize a list of (geometry, fill_value) tuples onto the given
	xarray coordinates. This only works for 1d latitude and longitude
	arrays.

	ARGUMENTS:
	----------
	shapes = [list] of tuples of (shapely.geom, fill_value)
	coords = [dict] of named 1d latitude and longitude arrays.
	latitude = [str] name of latitude key. default:'latitude'
	longitude = [str] name of longitude key. default:'longitude'
	fill = fill_value

	RETURNS:
	--------

	xarray.DataArray

	'''
	from rasterio import features
	import xarray as xr
	if fill == None:
		fill = np.nan
	transform = transform_from_latlon( coords[ latitude ], coords[ longitude ] )
	out_shape = ( len( coords[ latitude ] ), len( coords[ longitude ] ) )
	raster = features.rasterize( shapes, out_shape=out_shape,
								fill=fill, transform=transform,
								dtype=float, **kwargs )
	# spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
	# return xr.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))
	return raster

def make_mask( fn ):
	with rasterio.open(fn) as rst:
		meta = rst.meta.copy()
		arr = np.empty_like(rst.read(1))

	# build the shape data we need in EPSG:3338
	shp = make_gdf()
	pts = shp.geometry.tolist()
	pts = [ (i,count+1) for count,i in enumerate(pts) ]
	return rasterize( pts, fill=0, out=arr, transform=meta['transform'], all_touched=True, dtype='float32' )

def open_raster( fn, band=1 ):
	with rasterio.open(fn) as rst:
		arr = rst.read(band)
	return arr

def extract_values( files, mask ):
	pool = mp.Pool(64)
	f = partial(open_raster, band=1)
	arr = np.array(pool.map( f, files ))
	pool.close()
	pool.join()
	pool = None; del pool
	mask_vals = np.unique(mask[mask > 0])
	out = dict()
	for mask_val in mask_vals:
		ind = zip(*np.where(mask == mask_val))
		for i,j in ind:
			out.update({ mask_val:arr[:,i,j] })
	del arr
	return out

if __name__ == '__main__':
	import os, glob, rasterio, itertools
	import pandas as pd
	import numpy as np
	import rasterio
	from rasterio.features import rasterize
	from shapely.geometry import Point
	import geopandas as gpd
	from functools import partial
	import multiprocessing as mp

	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	output_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/project_data_delivery/Keith_DOT_extractions'
	template_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/5ModelAvg/historical/tasmin/tasmin_mean_C_ar5_5ModelAvg_historical_08_2004.tif'
	models = ['CRU-TS40']
	scenarios = ['historical']
	variables = ['pr']
	
	# make mask
	mask = make_mask( template_fn )
	output_dict = dict()
	for model, scenario, variable in itertools.product(models, scenarios, variables):
		cur_path = os.path.join(base_dir,model,scenario,variable)
		files_df = list_data( cur_path ) # these are sorted
		decade_grouper = files_df.apply(lambda x:str(x.year)[:3], axis=1)
		file_groups = [ j.fn.tolist() for i,j in files_df.groupby( decade_grouper ) ]
		out = [ extract_values( files, mask ) for files in file_groups ]
		out_df = pd.concat([ pd.DataFrame(i) for i in out ]) # stack the file groups chronologically
		out_key = '{}_{}_{}'.format( model, scenario, variable )
		output_dict[ out_key ] = out_df
		print( 'completed:{}'.format(out_key) )

# # future index
# future_dates = pd.date_range('2006-01','2101-01',freq='M')
# future_dates = [ [str(i.year), str(i.month)] for i in future_dates ]
# future_dates = [ '-'.join([y,'0'+m]) if len(m) == 1 else '-'.join([y,m]) for y,m in future_dates  ]

# historical index -- data needs slicing...
historical_dates = pd.date_range('1901-01','2016-01',freq='M')
historical_dates = [[str(i.month), str(i.year)] for i in historical_dates ]
historical_dates = [ '-'.join([y,'0'+m]) if len(m) == 1 else '-'.join([y,m]) for y,m in historical_dates ]

# make data frames
# df1_future, df2_future = [pd.DataFrame({key:np.array(output_dict[key][i]) for key in output_dict if 'historical' not in key }, index=future_dates) for i in [1,2]] 
df1_historical, df2_historical = [pd.DataFrame({key:np.array(output_dict[key][i])[-len(historical_dates):] for key in output_dict if 'historical' in key }, index=historical_dates) for i in [1,2]]

# dump them to disk
naming_lookup = {1:'LongLake', 2:'GlitterGulch'}

# df1_future_fn = 'precipitation_cmip5_allmodels_allscenarios_futures_2006-2100_LongLake_AK.csv'
# df2_future_fn = 'precipitation_cmip5_allmodels_allscenarios_futures_2006-2100_GlitterGulch_AK.csv'
df1_historical_fn = 'precipitation_cru_ts40_allmodels_allscenarios_historical_1901-2015_LongLake_AK.csv'
df2_historical_fn = 'precipitation_cru_ts40_allmodels_allscenarios_historical_1901-2015_GlitterGulch_AK.csv'

# df1_future.to_csv( os.path.join( output_dir, df1_future_fn), sep=',' )
# df2_future.to_csv( os.path.join( output_dir, df2_future_fn), sep=',' )
df1_historical.to_csv( os.path.join( output_dir, df1_historical_fn), sep=',' )
df2_historical.to_csv( os.path.join( output_dir, df2_historical_fn), sep=',' )

