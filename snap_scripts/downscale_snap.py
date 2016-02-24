# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
import glob, os
import downscale

# SETUP BASELINE
clim_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = downscale.Baseline( filelist )

# SETUP DATASET
output_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/test'
future_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/snap_prepped_data/IPSL-CM5A-LR/hur/hur_Amon_IPSL-CM5A-LR_rcp26_r1i1p1_200601_210012.nc'
historical_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/snap_prepped_data/IPSL-CM5A-LR/hur/hur_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001_200512.nc'
variable = 'hur'
model = 'IPSL-CM5A-LR'
scenario = 'rcp26'
historical = downscale.Dataset( historical_fn, variable, model, scenario, units=None )
future = downscale.Dataset( future_fn, variable, model, scenario, units=None )

# DOWNSCALE
clim_begin = '1961'
clim_end = '1990'
ar5 = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, future, \
		metric='mean', ds_type='absolute', level=1000, level_name='plev' )
ar5.downscale( output_dir=output_dir )

# CRU historical
output_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/test'
historical_fn = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.cld.dat.nc'
historical = downscale.Dataset( historical_fn, 'cld', 'cru_ts31', 'observed', units=None, interp=True )
clim_begin = '1961'
clim_end = '1990'
cru = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, metric='mean', ds_type='relative' )
cru.downscale( output_dir=output_dir )

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 