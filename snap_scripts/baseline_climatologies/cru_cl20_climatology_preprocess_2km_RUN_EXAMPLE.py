# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # EXAMPLE RUN OF THE ABOVE FOR TEM DATA CREATION
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# REGRID CRU-CL20 to SNAP AKCAN 2km
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
import subprocess, os, glob

script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/baseline_climatologies'
os.chdir( script_path )
base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
cru_filenames = glob.glob('/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/*.dat.gz')
template_raster_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/templates/akcan_2km/tas_mean_C_ar5_IPSL-CM5A-LR_rcp26_01_2006.tif'
for cru_filename in cru_filenames:
	print( 'working on: {}'.format( os.path.basename( cru_filename ) ) )
	variable = os.path.basename( cru_filename ).split( '_' )[-1].split('.')[0]
	done = subprocess.call([ 'ipython', 'cru_cl20_climatology_preprocess_2km.py', '--','-p', base_path, '-cru', cru_filename, '-v', variable ,'-tr', template_raster_fn ])

