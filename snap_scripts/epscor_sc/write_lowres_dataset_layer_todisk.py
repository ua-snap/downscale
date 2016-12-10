# # WRITE TO DISK
# tmp = ar5.ds.to_dataset('pre')
# tmp.to_netcdf( out_fn.replace('.tif', '.nc') )

import xarray as xr
import numpy as np
import rasterio, os

class WriteDataset( object ):
	def __init__( self, fn, variable, lat='lat', lon='lon', northup=False ):
		self.fn = fn
		self.variable = variable
		self.lat = lat
		self.lon = lon
		self.northup = northup
		self.ds = self.ds()

		if self.northup == True:
			self._northup( latitude = self.lat )

	def ds( self ):
		ds = xr.open_dataset( self.fn )
		return ds
	def _northup( self, latitude='lat' ):
		''' this works only for global grids to be downscaled flips it northup '''
		if self.ds[ latitude ][0].data < 0: # meaning that south is north globally
			self.ds[ latitude ] = np.flipud( self.ds[ latitude ] )
			# flip each slice of the array and make a new one
			flipped = np.array( [ np.flipud( arr ) for arr in self.ds[ self.variable ].data ] )
			self.ds[ self.variable ] = (('time', 'lat', 'lon' ), flipped )
	@property
	def affine( self ):
		''' simple way to make an affine transform from lats and lons coords '''
		from affine import Affine
		lats = np.asarray( self.ds[ self.lat ] )
		lons = np.asarray( self.ds[ self.lon ] )
		if (np.max( lats ) - 90) < np.abs( np.mean( np.diff( lats ) ) ):
			lat_max = 90.0
		else:
			lat_max = np.max( lats )

		# set the lonmax to the corner. --> this can get you into trouble with non-global data
		# but I am unsure how to make it more dynamic at the moment. [ML]
		lon_arr = np.array([ -180.0, 0.0 ])
		idx = (np.abs(lon_arr - np.min( lons ) ) ).argmin()
		lon_max = lon_arr[ idx ]

		trans = Affine.translation(lon_max, lat_max)
		scale = Affine.scale(lons[1] - lons[0], lats[1] - lats[0])
		return trans * scale		
	def to_gtiff( self, output_filename, layers=[0] ):
		''' write the mask to geotiff given an output_filename '''
		meta = {'compress':'lzw'}
		arr = np.array([ self.ds[self.variable][layer].data for layer in layers ])
		count, height, width = arr.shape
		meta.update( affine=self.affine, height=height, width=width, 
					count=arr.shape[0], dtype='float32', driver='GTiff' )
		with rasterio.open( output_filename, 'w', **meta ) as out:
			out.write( arr.astype( np.float32 ) )
