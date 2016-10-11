# THIS IS MY EXPLANATION OF HOW I HAVE ACTUALLY PRODUCED THE CORRECT DATA 


# FIRST 
'''
I used the following code to generate a single output raster for January 2006 anomaly created from files
we have already found to be identical to the anomalies Matt has created in the past. 

This is a test of `R` package raster projectRaster function, which wraps some of gdalwarp

'''

```R
setwd( "/workspace/Shared/Users/mfleonawicz/tmpDir" )
library(raster)

# Michael's data
b.clim <- brick("/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/matt_mike_test/tas_GFDL-CM3_rcp60_multiband_climatology.tif")
b.anom <- brick("/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/matt_mike_test/tas_GFDL-CM3_rcp60_multiband_anomalies.tif")
trial2km <- readAll(raster("/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_minmax/GFDL-CM3/rcp60/tas/anom/tas_mean_C_ar5_GFDL-CM3_rcp60_01_2041_anom.tif"))
yrs <- rep(2006:2100, each=12)
anom0 <- flip(rotate(subset(b.anom, which(yrs==2006))), "y")
clim0 <- flip(rotate(b.clim), "y")
r <- subset( anom0, 1 )
ext <- extent(-180, 180, -90, 90)
r@extent <- ext
test2km <- projectRaster(r, trial2km)
writeRaster(test2km, "TMP_test2km_2006-01_snapped.tif", datatype="FLT4S", overwrite=TRUE)
```


# SECOND 
'''
I ran a GDALWARP version using the python code below and sending functions to system commands
'''

```python

import os, rasterio
import numpy as np

# # generate a raster to project into using the gdalwarp CLI
fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled/GFDL-CM3/rcp60/tas/tas_mean_C_ar5_GFDL-CM3_rcp60_01_2006.tif'
with rasterio.open( fn ) as out:
	meta = out.meta
	arr = out.read( 1 )

meta.update( compress='lzw' )
meta.pop( 'transform' ) # drop the old-style geotransform

output_filename = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/testwarp/tas_mean_C_ar5_GFDL-CM3_rcp60_01_2006_anom.tif' 
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write( np.zeros_like( arr ), 1 )

# # # # # 

# write out the first slice of the same anomalies as above to use in test
anom = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/matt_mike_test/tas_GFDL-CM3_rcp60_multiband_anomalies.tif'
with rasterio.open( anom ) as rst:
	meta = rst.meta
	arr = rst.read( 1 ) # band 1 here is Jan 2006 

meta.update( compress='lzw', count=1 )
meta.pop( 'transform' )

anom_out = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/testwarp/tas_mean_C_ar5_GFDL-CM3_rcp60_01_2006_anom_LR.tif'
with rasterio.open( anom_out, 'w', **meta ) as rst:
	rst.write( arr, 1 )

anom_out_snap = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/testwarp/tas_mean_C_ar5_GFDL-CM3_rcp60_01_2006_anom_LR_SNAPPED.tif'
anom_out_snap_new = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/testwarp/tas_mean_C_ar5_GFDL-CM3_rcp60_01_2006_anom_LR_SNAPPED_FLIPPED.tif'

# gdal_translate command to reset the bounds to the proper upper right corners in the acceptable space.
#  note: these data are south-up and is shown as a_ullr below
done = os.system( ' '.join([ 'gdal_translate', '-a_ullr','0 -90 360 90', anom_out, anom_out_snap ]) )
# warp them to north up -- I have checked this and it works correctly
done = os.system( ' '.join([ 'gdalwarp', '-overwrite', '-te','0 90 360 -90', anom_out_snap, anom_out_snap_new ]) )
# warp into the empty raster we made above using a BILINEAR interpolation
done = os.system( ' '.join([ 'gdalwarp', '-r', 'BILINEAR', anom_out_snap_new, output_filename ]) )

```

# THIRD

'''
I grabbed the regridding code from the package and ran the same data through this using the python code below:

'''

```python

def interp_ds( anom, base, src_crs, src_nodata, dst_nodata, src_transform, resample_type='bilinear',*args, **kwargs ):
	'''	
	anom = [numpy.ndarray] 2-d array representing a single monthly timestep of the data to be downscaled. 
							Must also be representative of anomalies.
	base = [str] filename of the corresponding baseline monthly file to use as template and downscale 
							baseline for combining with anomalies.
	src_transform = [affine.affine] 6 element affine transform of the input anomalies. [should be greenwich-centered]
	resample_type = [str] one of ['bilinear', 'count', 'nearest', 'mode', 'cubic', 'index', 'average', 'lanczos', 'cubic_spline']
	'''	
	import rasterio
	from rasterio.warp import reproject, RESAMPLING

	resampling = {'average':RESAMPLING.average,
				'cubic':RESAMPLING.cubic,
				'lanczos':RESAMPLING.lanczos,
				'bilinear':RESAMPLING.bilinear,
				'cubic_spline':RESAMPLING.cubic_spline,
				'mode':RESAMPLING.mode,
				'count':RESAMPLING.count,
				'index':RESAMPLING.index,
				'nearest':RESAMPLING.nearest }
	
	base = rasterio.open( base )
	baseline_arr = base.read( 1 )
	baseline_meta = base.meta
	baseline_meta.update( compress='lzw' )
	output_arr = np.zeros_like( baseline_arr ) # THIS COULD BE DANGEROUS?
	reproject( anom, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
			dst_transform=baseline_meta['affine'], dst_crs=baseline_meta['crs'],\
			dst_nodata=dst_nodata, resampling=resampling[ resample_type ], SOURCE_EXTRA=1000 )
	return output_arr

if __name__ == '__main__':
	import rasterio
	import numpy as np

	# use the file we snapped and flipped globally above for the GDALWARP example
	fn = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/testwarp/tas_mean_C_ar5_GFDL-CM3_rcp60_01_2006_anom_LR_SNAPPED_FLIPPED.tif'
	base = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/testwarp/tas_mean_C_ar5_GFDL-CM3_rcp60_01_2006_anom.tif'
	meta = rasterio.open( base ).meta 
 	meta.update( compress='lzw' )

	with rasterio.open( fn ) as anom:
		anom_meta = anom.meta
		anom = anom.read( 1 )
		src_crs = anom_meta[ 'crs' ]
		src_nodata = None
		dst_nodata = None
		src_transform = anom_meta[ 'affine' ]

	new_arr = interp_ds( anom, base, src_crs, src_nodata, dst_nodata, src_transform )

	with rasterio.open( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/testwarp/warped_3338_from_PCLL.tif', 'w', **meta ) as out:
		out.write( new_arr, 1 )

```

# CONCLUSIONS

'''
What I have found is that the data created using projectRaster and the data created using rasterio.warp are identical, with the small 
exception of the dateline in the projectRaster version that can be overcome, but doesnt affect the land-based data we provide. The same 
differences exist between the commandline gdalwarp and both projectRaster and rasterio.warp versions. This means that whatever differences we
are seeing between the old and new tas versions is attributable to something else, if indeed the regridding was performed using projectRaster
and not something else. So then I thought, lets just take these files I have just downscaled and remove the anomaly from that downscaled January
2006 file to make SURE that what we have in return is PRISM.  This is INDEED the case! and is perfectly generated.  Now I am REALLY confused, 
as it seems that what I have done is 100% correct. Which would either mean I am missing something completely, or there is something wrong with 
the old downscaled tas data...  I realize this means more investigating, but we need to figure out why this is happening, since it appears that
the downscale package is doing exactly what it is supposed to, and is the same as its brother from another mother projectRaster / R.  
'''