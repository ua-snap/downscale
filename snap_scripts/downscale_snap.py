# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# BASELINE:
import glob, os
import downscale
clim_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = downscale.Baseline( filelist )

# DATASET:
fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/snap_prepped_data/IPSL-CM5A-LR/hur/hur_Amon_IPSL-CM5A-LR_rcp26_r1i1p1_200601_210012.nc'
variable = 'hur'
model = 'IPSL-CM5A-LR'
scenario = 'rcp26'

# ar5 modeled
ds = Dataset( fn, variable, model, scenario, units=None )
clim_begin = '2061'
clim_end = '2090'
ar5 = downscale.DeltaDownscale( baseline, clim_begin, clim_end, ds, future=None, \
		metric='mean', ds_type='absolute', level=1000, level_name='plev' )
ar5.downscale( output_dir='/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/test' )

# CRU historical
cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.cld.dat.nc'
cru_ds = downscale.Dataset( cru_ts, 'cld', 'cru_ts31', 'observed', units=None, interp=True )
clim_begin = '1961'
clim_end = '1990'
cru = DeltaDownscale( baseline, clim_begin, clim_end, cru_ds, metric='mean', ds_type='relative' )
cru.downscale( output_dir='/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/test' )

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 