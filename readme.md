[![Stories in Ready](https://badge.waffle.io/ua-snap/downscale.svg?label=ready&title=Ready)](http://waffle.io/ua-snap/downscale)

downscale
---------

Python Package for Delta-Downscaling of Climatic Research Unit's TS 3.x OR PCMDI's CMIP5 data to a set of baseline 
monthly climatologies. In (potentially) different (overlapping) reference system, resolution, area of interest.

#### Note:
This package is not yet complete and is expected NOT to work as expected until full release.  It is also geared specifically
for the needs of Scenarios Network for Alaska + Arctic Planning (SNAP).

#### EXAMPLES

##### AR5
```python
if __name__ == '__main__':
	from downscale import DownscaleAR5
	# input args
	ar5_modeled = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped/IPSL-CM5A-LR/clt/clt_Amon_IPSL-CM5A-LR_rcp26_r1i1p1_200601_210012.nc'
	ar5_historical = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped/IPSL-CM5A-LR/clt/clt_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001_200512.nc'
	clim_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
	template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	base_path = '/atlas_scratch/malindgren/CMIP5/TEST_AR5'

	# EXAMPLE RUN
	down = DownscaleAR5.DownscaleAR5( ar5_modeled, ar5_historical, base_path, clim_path, template_raster_fn=template_raster_fn, ncores=32 ) #, climatology_begin, climatology_end, plev, absolute, metric, ncores )
	output = down.downscale_ar5_ts()
```

##### CRU
```python
if __name__ == '__main__':
	from downscale import DownscaleCRU
	
	# example of post_downscale_function - pass in at DownscaleCRU()
	def clamp_vals( x ):
		''' clamp the values following the relative humidity downscaling '''
		x[ (x > 100) & (x < 500) ] = 95
		return x

	# minimum required arguments
	cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.cld.dat.nc'
	clim_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
	template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	base_path = '/atlas_scratch/malindgren/CMIP5/CRU2'
	
	# run example
	down = DownscaleCRU.DownscaleCRU( cru_ts, clim_path, template_raster_fn, base_path, absolute=False, ncores=32 )
	output = down.downscale_cru_ts()
```