# CRU TS3.1 Run -- January 2016

import glob, os, itertools, rasterio
from downscale import DeltaDownscale, Baseline, Dataset

# static args setup
cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS31/cru_ts_3_10.1901.2009.tmx.dat.nc.gz'
clim_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2/tmax'
output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru/CRU_TS31/tmx' #'/Data/malindgren/downscale_outputs/CRU'
ncores = 32
clim_begin = '1961'
clim_end = '1990'
variable = 'tmx'
model = 'cru_ts31'
scenario = 'historical'

# RUN 2.0
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = Baseline( filelist )

# DOWNSCALE
mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
clim_begin = '1961'
clim_end = '1990'

# FOR CRU WE PASS THE interp=True so we interpolate across space first when creating the Dataset()
historical = Dataset( cru_ts, variable, model, scenario, units=None, interp=True, ncpus=32 )

ar5 = DeltaDownscale( baseline, clim_begin, clim_end, historical, future=None, \
		metric='mean', downscaling_operation='add', mask=mask, mask_value=0, ncpus=32, \
		src_crs={'init':'epsg:4326'}, src_nodata=-9999, dst_nodata=None,
		post_downscale_function=None )

ar5.downscale( output_dir=output_path )


# OLDER STUFF
# addition for later?
# if variable == 'hur':
# 	def clamp_vals( x ):
# 		''' clamp the values following the relative humidity downscaling '''
# 		x[ (x > 100) & (x < 500) ] = 95
# 		return x