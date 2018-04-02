# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# MAKE GREENWICH CENTERED AOI_MASK FOR CRU 10min RUNS THAT WE USE
# WHEN CORRECTING PRECIP VALUES THAT ARE WAAAY TOO HIGH.
# 
# Michael Lindgren (April 2018)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

def shiftgrid( lon0, datain, lonsin, start=True, cyclic=360.0 ):
	"""
	Shift global lat/lon grid east or west.
	.. tabularcolumns:: |l|L|
	==============   ====================================================
	Arguments        Description
	==============   ====================================================
	lon0             starting longitude for shifted grid
					 (ending longitude if start=False). lon0 must be on
					 input grid (within the range of lonsin).
	datain           original data with longitude the right-most
					 dimension.
	lonsin           original longitudes.
	==============   ====================================================
	.. tabularcolumns:: |l|L|
	==============   ====================================================
	Keywords         Description
	==============   ====================================================
	start            if True, lon0 represents the starting longitude
					 of the new grid. if False, lon0 is the ending
					 longitude. Default True.
	cyclic           width of periodic domain (default 360)
	==============   ====================================================
	returns ``dataout,lonsout`` (data and longitudes on shifted grid).
	"""
	if np.fabs(lonsin[-1]-lonsin[0]-cyclic) > 1.e-4:
		# Use all data instead of raise ValueError, 'cyclic point not included'
		start_idx = 0
	else:
		# If cyclic, remove the duplicate point
		start_idx = 1
	if lon0 < lonsin[0] or lon0 > lonsin[-1]:
		raise ValueError('lon0 outside of range of lonsin')
	i0 = np.argmin(np.fabs(lonsin-lon0))
	i0_shift = len(lonsin)-i0 # [ML] THIS COULD BE THE LINE GETTING US!!!
	if np.ma.isMA(datain):
		dataout  = np.ma.zeros(datain.shape,datain.dtype)
	else:
		dataout  = np.zeros(datain.shape,datain.dtype)
	if np.ma.isMA(lonsin):
		lonsout = np.ma.zeros(lonsin.shape,lonsin.dtype)
	else:
		lonsout = np.zeros(lonsin.shape,lonsin.dtype)
	if start:
		lonsout[0:i0_shift] = lonsin[i0:]
	else:
		lonsout[0:i0_shift] = lonsin[i0:]-cyclic
	dataout[...,0:i0_shift] = datain[...,i0:]
	if start:
		lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]+cyclic
	else:
		lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]
	dataout[...,i0_shift:] = datain[...,start_idx:i0+start_idx]
	return dataout,lonsout
def rotate( dat, lons, to_pacific=False ):
	'''rotate longitudes in WGS84 Global Extent'''
	if to_pacific == True:
		# to 0 - 360
		dat, lons = shiftgrid( 0., dat, lons )
	elif to_pacific == False:
		# to -180.0 - 180.0 
		dat, lons = shiftgrid( 180., dat, lons, start=False )
	else:
		raise AttributeError( 'to_pacific must be boolean True:False' )
	return dat, lons

def transform_from_latlon( lat, lon ):
	''' simple way to make an affine transform from lats and lons coords '''
	from affine import Affine
	lat = np.asarray( lat )
	lon = np.asarray( lon )
	trans = Affine.translation(lon[0], lat[0])
	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
	return trans * scale
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
def bounds_to_extent( bounds ):
	'''
	take input rasterio bounds object and return an extent
	'''
	l,b,r,t = bounds
	return [ (l,b), (r,b), (r,t), (l,t), (l,b) ]

if __name__ == '__main__':
	import os
	import xarray as xr
	import geopandas as gpd
	import numpy as np
	import rasterio
	from shapely.geometry import Polygon

	# get raw cru ts40
	fn = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS40/cru_ts4.00.1901.2015.pre.dat.nc.gz'
	shp_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/masks/pcll_template_10min_extent_with_nwt.shp'

	# open it
	ds = xr.open_dataset(fn)
	shp = gpd.read_file( shp_fn )

	# rotate it to pcll
	dat = np.flipud(ds.pre[0,...].data)
	lons = np.array(ds.lon)
	
	dat_pc, lons_pc = rotate( dat, lons, to_pacific=True )

	coords = {'lon':lons_pc, 'lat':np.flipud(np.array(ds.lat))}

	# rasterize using the pcll file with 
	shapes = [ (i,1) for i in shp.geometry ]
	rst = rasterize( shapes, coords, latitude='lat', longitude='lon', fill=0 )

	# rotate back to GCLL
	dat_gc, lons_gc = rotate( rst, lons_pc, to_pacific=False )
	height,width = dat_gc.shape
	meta = {
		'driver':'GTiff',
		'height':height,
		'width':width,
		'count':1,
		'crs':{'init':'epsg:4326'},
		'dtype':'float32',
		'transform':transform_from_latlon( coords[ 'lat' ], lons_gc ),
		'compress':'lzw'
	}

	with rasterio.open( '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/TEST_GCLL_CRU.tif', 'w', **meta ) as tmp:
		tmp.write( dat_gc.astype(np.float32), 1 )

	# polygonize and make a polygon from the bounds...
	new_ext,val = [ i for i in rasterio.features.shapes( dat_gc.astype(np.float32), mask=dat_gc==1, transform=meta['transform']) ][0]
	pol = Polygon( bounds_to_extent( Polygon(new_ext['coordinates'][0]).bounds ) )

	# make a geodataframe and dump out as a shapefile
	new_df = gpd.GeoDataFrame({'id':[1],'geometry':[pol]}, crs={'init':'epsg:4326'}, geometry='geometry' )
	new_df.to_file( shp_fn.replace('pcll','gcll') )

