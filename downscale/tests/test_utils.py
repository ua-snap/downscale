# -*- coding: utf8 -*-
# # # # 
# beginnings of a test suite for the utility functions
# # # # 

import unittest
import rasterio

class TestWriteGTiff( unittest.TestCase ):
	''' tests for testing write_gtiff function '''
	def setUp( self ):
		print( 'setup' )
		import rasterio
		from affine import Affine
		import numpy as np

		# make dummy rasters / arrays for testing
		meta = {'affine': Affine(1.0, 0.0, -180.0, 0.0, -1.0, 180.0),
				'count': 1,
				'crs': {'init': 'epsg:4326'},
				'driver': u'MEM',
				'dtype': 'float32',
				'height': 180,
				'nodata': -3.4e+38,
				'width': 360}

		with rasterio.drivers():
			self.rst = rasterio.open( '', mode='w', **meta )
			self.arr = self.rst.read( 1 )
			meta.update( driver=u'GTiff')
			self.meta = meta
			self.test_fn = './test_write_gtiff_59137.tif'
	def test_output_crs( self ):
		print( 'test_meta' )
		import rasterio
		from downscale.DownscalingUtils import write_gtiff
		out = write_gtiff( self.arr, self.meta, self.test_fn )
		new_meta = rasterio.open( out ).meta
		self.assertEqual( self.meta['crs'], new_meta['crs'] )
	def test_output_shape( self ):
		print( 'test_shape' )
		import rasterio
		from downscale.DownscalingUtils import write_gtiff
		out = write_gtiff( self.arr, self.meta, self.test_fn )
		new_meta = rasterio.open( out ).meta
		self.assertEqual( self.meta['height'], new_meta['height'] )
		self.assertEqual( self.meta['width'], new_meta['width'] )
	def tearDown( self ):
		print( 'teardown' )
		import os
		if os.path.exists( self.test_fn ):
			os.unlink( self.test_fn )

if __name__ == '__main__':
	unittest.main()
