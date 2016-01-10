from downscale import DownscalingUtils
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

rst = rasterio.open( '', mode='w', **meta )
arr = rst.read( 1 )

# test functions all options



