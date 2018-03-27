# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # EXAMPLE RUN OF THE ABOVE FOR TEM DATA CREATION
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # REGRID CRU-CL20 to SNAP AKCAN 10min
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
import subprocess, os, glob

script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/baseline_climatologies'
os.chdir( script_path )
base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
cru_filenames = glob.glob('/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/*.dat.gz')
template_raster_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/templates/akcan_10min/akcan_with_nwt_15k_template.tif'
for cru_filename in cru_filenames:
	print( 'working on: {}'.format( os.path.basename( cru_filename ) ) )
	variable = os.path.basename( cru_filename ).split( '_' )[-1].split('.')[0]
	done = subprocess.call([ 'ipython', 'cru_cl20_climatology_preprocess_10min.py', '--','-p', base_path, '-cru', cru_filename, '-v', variable ,'-tr', template_raster_fn ])
