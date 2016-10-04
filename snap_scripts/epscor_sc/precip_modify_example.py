# NEW PR APPROACH FOR CRU / CMIP5
def transform_from_latlon( lat, lon ):
	''' simple way to make an affine transform from lats and lons coords '''
	from affine import Affine
	lat = np.asarray( lat )
	lon = np.asarray(lon)
	trans = Affine.translation(lon[0], lat[0])
	scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
	return trans * scale

def find_boundary( arr ):
	'''
	return a mask of the boundary limit of DATA cells (overlays edge DATA not NA)
	this is especially useful if the data are land-only.  As in the CRU TS3.x data.
	'''
	bool_arr = np.copy( arr )
	ind = np.isnan( bool_arr )
	bool_arr[ ~ind ] = 1
	bool_arr[ ind ] = 0
	return find_boundaries( bool_arr, mode='inner' )

def correct_boundary( arr, bound_mask ):
	''' correct the boundary pixels with non-acceptable values '''
	upperthresh = np.percentile( arr[~np.isnan( arr )], 95 )
	ind = np.where( bound_mask == True )
	vals = arr[ ind ]
	vals[ vals < 0.5 ] = 0.5
	vals[ vals > upperthresh ] = upperthresh
	arr[ ind ] = vals
	return arr

def correct_inner( arr, bound_mask ):
	''' correct the inner pixels with non-acceptable values '''
	upperthresh = np.percentile( arr[~np.isnan( arr )], 95 )
	mask = np.copy( arr )	
	ind = np.where( (arr > 0) & bound_mask != True )
	vals = arr[ ind ]
	vals[ vals < 0.5 ] = np.nan # set to the out-of-bounds value
	vals[ vals > upperthresh ] = upperthresh 
	arr[ ind ] = vals
	return arr

def run( arr, bound_mask ):		
	arr = correct_boundary( arr, bound_mask )
	return correct_inner( arr, bound_mask )

if __name__ == '__main__':
	import xarray as xr
	import numpy as np
	from skimage.segmentation import find_boundaries
	from affine import Affine
	import rasterio

	# fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/cmip5/prepped/GFDL-CM3/historical/pr/pr_GFDL-CM3_historical_r1i1p1_1860_2005.nc'
	fn = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.pre.dat.nc'

	varname = 'pre' # 'pr'
	convert_to_mm = False
	landonly = True

	# read it
	ds = xr.open_dataset( fn )

	# calculate climatologies
	climatology = ds[ varname ].sel( time=slice('1961', '1990') ) 

	# CONVERT THE PRECIP TO MM/MONTH
	if convert_to_mm:
		days_in_month = 31
		climatology = climatology * 86400 * days_in_month

	# make the monthly means	
	climatology = climatology.groupby( 'time.month' ).mean( 'time' )

	# this is not working as it should since we need a way to deal with
	# CRU land-only data and the CMIP5 global grids...
	if landonly:
		# find the boundary -- where there *is* data
		bound_mask = find_boundary( climatology[ 1, ... ].data )

	# fix the data
	clim_fixed = np.array( [ np.flipud( run( i.data, bound_mask ) ) for i in climatology ] )
	
	# # WRITE SOME TEST FILES TO DISK FOR VISUAL ASSESSMENT
	affine = transform_from_latlon( np.flipud( ds.lat ), ds.lon )
	count, height, width = climatology.shape

	meta = { 'affine' : affine,
			'count' : count,
			'height' : height,
			'width' : width,
			'dtype' : 'float32',
			'driver' : 'GTiff' }

	## BOUNDARY FILE
	# output_filename = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/pr_test/bounday_test_inner.tif'
	# meta.update( count=1 )
	# with rasterio.open( output_filename,'w', **meta ) as out:
	# 	out.write( np.flipud(bound_mask).astype(np.float32), 1 )

	## MODIFIED DATA
	output_filename = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/pr_test/dat_test.tif'
	with rasterio.open( output_filename,'w', **meta ) as out:
		clim_fixed[ np.isnan(clim_fixed) ] = -99
		out.write( clim_fixed.astype( np.float32 ) )

