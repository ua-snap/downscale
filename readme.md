[![Stories in Ready](https://badge.waffle.io/ua-snap/downscale.svg?label=ready&title=Ready)](http://waffle.io/ua-snap/downscale)
downscale
---------

start of a readme for this new package

This package is not yet complete and is expected NOT to work as expected until full release

### EXAMPLES

##### AR5
	if __name__ == '__main__':
		# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		# example of use of the new DownscaleAR5 / DownscalingUtils classes
		# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		from downscale import DownscaleAR5
		# input args
		ar5_modeled = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped/IPSL-CM5A-LR/clt/clt_Amon_IPSL-CM5A-LR_rcp26_r1i1p1_200601_210012.nc'
		ar5_historical = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped/IPSL-CM5A-LR/clt/clt_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001_200512.nc'
		clim_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
		template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
		base_path = '/atlas_scratch/malindgren/CMIP5/TEST_AR5'

		# EXAMPLE RUN -- TESTING
		down = DownscaleAR5.DownscaleAR5( ar5_modeled, ar5_historical, base_path, clim_path, template_raster_fn=template_raster_fn, ncores=32 ) #, climatology_begin, climatology_end, plev, absolute, metric, ncores )
		output = down.downscale_ar5_ts()