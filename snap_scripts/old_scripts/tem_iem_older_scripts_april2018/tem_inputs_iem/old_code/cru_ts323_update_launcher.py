
# SCRIPT TO RUN THE CRU TS3.1 BUILT SCRIPT OVER THE CRU TS3.2.3 UPDATE (SEPT.2015)
# WHICH EXTENDS THE SERIES TO 12/2014.
# THIS IS CURRENTLY WORKING FOR CLD, TMP, VAP, AND MORE TO COME!
# # # # #
# Author: Michael Lindgren (malindgren@alaska.edu)
# # # # # 
# CURRENTLY SETUP TO RUN ON EOS.

# CLD
import os
os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
ncores = '14'
base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323'
cru_ts31 = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323/cru_ts3.23.1901.2014.cld.dat.nc'
cl20_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
anomalies_calc_type = 'relative'
downscaling_operation = 'mult'

climatology_begin = '1961'
climatology_end = '1990'
year_begin = '1901'
year_end = '2014'
variable = 'cld'
metric = 'pct'

args_tuples = [ ('hi', cru_ts31), ('ci', cl20_path), ('tr', template_raster_fn), 
				('base', base_path), ('bt', year_begin), ('et', year_end),
				('cbt', climatology_begin), ('cet', climatology_end),
				('nc', ncores), ('at', anomalies_calc_type), ('m', metric), 
				('dso', downscaling_operation), ('v', variable) ]

args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
os.system( 'ipython2.7 -- tas_cld_cru_ts31_to_cl20_downscaling.py ' + args )


# TAS 
import os
os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
ncores = '14'
base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323'
cru_ts31 = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323/cru_ts3.23.1901.2014.tmp.dat.nc'
cl20_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
anomalies_calc_type = 'absolute'
downscaling_operation = 'add'

climatology_begin = '1961'
climatology_end = '1990'
year_begin = '1901'
year_end = '2014'
variable = 'tas'
metric = 'C'

args_tuples = [ ('hi', cru_ts31), ('ci', cl20_path), ('tr', template_raster_fn), 
				('base', base_path), ('bt', year_begin), ('et', year_end),
				('cbt', climatology_begin), ('cet', climatology_end),
				('nc', ncores), ('at', anomalies_calc_type), ('m', metric), 
				('dso', downscaling_operation), ('v', variable) ]

args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
os.system( 'ipython2.7 -- tas_cld_cru_ts31_to_cl20_downscaling.py ' + args )


# VAP (HUR)
import os
os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
ncores = '14'
base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323'
cru_ts31_vap = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323/cru_ts3.23.1901.2014.vap.dat.nc'
cru_ts31_tas = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323/cru_ts3.23.1901.2014.tmp.dat.nc'
cl20_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/hur/akcan' # hur 
template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
anomalies_calc_type = 'relative'
downscaling_operation = 'mult'

climatology_begin = '1961'
climatology_end = '1990'
year_begin = '1901'
year_end = '2014'
variable = 'hur'
metric = 'pct'

args_tuples = [ ('hhi', cru_ts31_vap), ('thi', cru_ts31_tas), ('ci', cl20_path), ('tr', template_raster_fn), 
				('base', base_path), ('bt', year_begin), ('et', year_end), 
				('cbt', climatology_begin), ('cet', climatology_end), 
				('nc', ncores), ('at', anomalies_calc_type), ('m', metric), 
				('dso', downscaling_operation), ('v', variable) ]

args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
os.system( 'ipython2.7 -i -- hur_cru_ts31_to_cl20_downscaling.py ' + args )

# PRECIP
import os
os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
ncores = '14'
base_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323'
cru_ts31 = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts323/cru_ts3.23.1901.2014.pre.dat.nc'
cl20_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/pre/akcan'
template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
anomalies_calc_type = 'relative'
downscaling_operation = 'mult'

climatology_begin = '1961'
climatology_end = '1990'
year_begin = '1901'
year_end = '2014'
variable = 'pre'
metric = 'mm'

args_tuples = [ ('hi', cru_ts31), ('ci', cl20_path), ('tr', template_raster_fn), 
				('base', base_path), ('bt', year_begin), ('et', year_end),
				('cbt', climatology_begin), ('cet', climatology_end),
				('nc', ncores), ('at', anomalies_calc_type), ('m', metric), 
				('dso', downscaling_operation), ('v', variable) ]

args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
os.system( 'ipython2.7 -- tas_cld_cru_ts31_to_cl20_downscaling.py ' + args )

