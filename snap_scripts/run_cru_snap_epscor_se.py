# CRU TS3.23 Run -- July 2016

# # # # TMX # # # # # # # 
import glob, os, itertools, rasterio
from downscale import DeltaDownscale, Baseline, Dataset, utils

# static args setup
cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.tmx.dat.nc'
clim_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/akcan/tmax'
output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_TEST/CRU_TS323/tmx'
ncores = 32
clim_begin = '1961'
clim_end = '1990'
variable = 'tmx'
model = 'ts323'
scenario = 'historical'
metric = 'mean'
out_varname = 'tasmax'
project = 'cru'
units = 'C'

# RUN 2.0
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
filelist = [ fn for fn in filelist if not '_14_' in fn ]
baseline = Baseline( filelist )

# DOWNSCALE
mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
clim_begin = '1961'
clim_end = '1990'

# FOR CRU WE PASS THE interp=True so we interpolate across space first when creating the Dataset()
historical = Dataset( cru_ts, variable, model, scenario, project, units, metric=metric, interp=True, method='linear', ncpus=32 )

ar5 = DeltaDownscale( baseline, clim_begin, clim_end, historical, future=None, \
		downscaling_operation='add', mask=mask, mask_value=0, ncpus=32, \
		src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
		post_downscale_function=None, varname=out_varname, modelname=None ) # -9999

if not os.path.exists( output_path ):
	os.makedirs( output_path )

ar5.downscale( output_dir=output_path )




# # # # OLD AND BROKEN! ! ! ! ! # # # # # # # 

# # # # TMN # # # # # # # 
import glob, os, itertools, rasterio
from downscale import DeltaDownscale, Baseline, Dataset, utils
		
# static args setup
cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.tmn.dat.nc'
clim_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism_v2/tmin'
output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru/CRU_TS323/tmn' #'/Data/malindgren/downscale_outputs/CRU'
ncores = 32
clim_begin = '1961'
clim_end = '1990'
variable = 'tmn'
model = 'ts323'
scenario = 'historical'
metric = 'mean'
out_varname = 'tasmin'
project = 'cru'
units = 'C'

# RUN 2.0
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = Baseline( filelist )

# DOWNSCALE
mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
clim_begin = '1961'
clim_end = '1990'

# FOR CRU WE PASS THE interp=True so we interpolate across space first when creating the Dataset()
historical = Dataset( cru_ts, variable, model, scenario, project=project, units=units, interp=True, ncpus=32 )

ar5 = DeltaDownscale( baseline, clim_begin, clim_end, historical, future=None, \
		metric=metric, downscaling_operation='add', mask=mask, mask_value=0, ncpus=32, \
		src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
		post_downscale_function=None, varname=out_varname, modelname=None ) #-9999

if not os.path.exists( output_path ):
	os.makedirs( output_path )

ar5.downscale( output_dir=output_path )

# # # # PRE # # # # # # # 
import glob, os, itertools, rasterio
from downscale import DeltaDownscale, Baseline, Dataset, utils

# static args setup
cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.pre.dat.nc'
clim_path = '/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_CAN_2km_PRISM/AK_CAN_geotiffs/pr/ak83albers'
output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru/CRU_TS323/pre' #'/Data/malindgren/downscale_outputs/CRU'
ncores = 32
clim_begin = '1961'
clim_end = '1990'
variable = 'pre'
model = 'ts323'
scenario = 'historical'
metric = 'total'
out_varname = 'pr'
project = 'cru'
units = 'mm'

# RUN 2.0
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = Baseline( filelist )

# DOWNSCALE
mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
clim_begin = '1961'
clim_end = '1990'

# FOR CRU WE PASS THE interp=True so we interpolate across space first when creating the Dataset()
historical = Dataset( cru_ts, variable, model, scenario, units, metric, interp=True, method='linear', ncpus=32 )

# new = interp_na( historical, 'cubic' )
ar5 = DeltaDownscale( baseline, clim_begin, clim_end, historical, future=None, \
		metric=metric, downscaling_operation='add', mask=mask, mask_value=0, ncpus=32, \
		src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
		post_downscale_function=None, varname=out_varname, modelname=None ) # -9999

if not os.path.exists( output_path ):
	os.makedirs( output_path )

ar5.downscale( output_dir=output_path )


# # # # TAS # # # # # # # 
import glob, os, itertools, rasterio
from downscale import DeltaDownscale, Baseline, Dataset, utils

# static args setup
cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.tmp.dat.nc'
clim_path = '/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_CAN_2km_PRISM/AK_CAN_geotiffs/tas/ak83albers'
output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru/CRU_TS323/tmp' #'/Data/malindgren/downscale_outputs/CRU'
ncores = 32
clim_begin = '1961'
clim_end = '1990'
variable = 'tmp'
model = 'ts323'
scenario = 'historical'
metric = 'mean'
out_varname = 'tas'
project = 'cru'
units = 'C'

# RUN 2.0
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = Baseline( filelist )

# DOWNSCALE
mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
clim_begin = '1961'
clim_end = '1990'

# FOR CRU WE PASS THE interp=True so we interpolate across space first when creating the Dataset()
historical = Dataset( cru_ts, variable, model, scenario, units, metric, interp=True, method='linear', ncpus=32 )

# new = interp_na( historical, 'cubic' )
ar5 = DeltaDownscale( baseline, clim_begin, clim_end, historical, future=None, \
		metric=metric, downscaling_operation='add', mask=mask, mask_value=0, ncpus=32, \
		src_crs={'init':'epsg:4326'}, src_nodata=None, dst_nodata=None,
		post_downscale_function=None, varname=out_varname, modelname=None ) # -9999

if not os.path.exists( output_path ):
	os.makedirs( output_path )

ar5.downscale( output_dir=output_path )





# OLDER STUFF
# addition for later?
# if variable == 'hur':
# 	def clamp_vals( x ):
# 		''' clamp the values following the relative humidity downscaling '''
# 		x[ (x > 100) & (x < 500) ] = 95
# 		return x