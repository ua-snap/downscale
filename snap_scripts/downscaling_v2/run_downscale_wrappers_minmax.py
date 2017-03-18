# run_downscale -- CMIP5 / CRU TS323 / CRU TS40
#  for variables: tas, pr, vap, hur, clt
#	*tasmin/tasmax require tas to be run first so we
#	perform the computation in a second run.
# # # # # # # # # # # # # # # # # # # # # # # # # # # 
import os, shutil

base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
# change the CCSM4 naming to NCAR-CCSM4 (if needed)
ccsm4_path = os.path.join( base_dir, 'downscaled', 'CCSM4' )
if os.path.exists( ccsm4_path ):
	_ = subprocess.call([ 'mv', path, path.replace( 'CCSM4', 'NCAR-CCSM4' ) ])
	print( 'changed CCSM4 error to NCAR-CCSM4 for proper min/max handling' )
del ccsm4_path	# end naming change

# change dir to the scripts directory
os.chdir( '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_v2' )

run_scripts = [ 'wrap_downscaler_cmip5_slurm_minmax.py',
				'wrap_downscaler_cru_40_slurm_minmax.py',
				'wrap_downscaler_cru_slurm_minmax.py']

for script in run_scripts:
	os.system( 'ipython {}'.format( script ) )

