# CRU TS3.1 Run -- January 2016

import glob, os, itertools
from downscale import DownscaleCRU

input_path = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323'

# static args setup
cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.cld.dat.nc'
clim_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
base_path = '/Data/malindgren/downscale_outputs/CRU'
ncores = 16

# list the files in the CRU NetCDF directory
nc_list = glob.glob( os.path.join( input_path, '*.nc' ) )

for nc in nc_list:
	print( nc )
	args = {}
	args.update( cru_ts=nc, clim_path=clim_path, template_raster_fn=template_raster_fn, base_path=base_path, ncores=ncores )

	# run
	down = DownscaleCRU.DownscaleCRU( nc, clim_path, template_raster_fn, base_path, ncores=ncores )
	output = down.downscale_cru_ts()



# addition for later?
# if variable == 'hur':
# 	def clamp_vals( x ):
# 		''' clamp the values following the relative humidity downscaling '''
# 		x[ (x > 100) & (x < 500) ] = 95
# 		return x