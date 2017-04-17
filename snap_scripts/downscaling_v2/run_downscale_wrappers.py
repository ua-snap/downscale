# run_downscale -- CMIP5 / CRU TS323 / CRU TS40
#  for variables: tas, pr, vap, hur, clt
#	*tasmin/tasmax require tas to be run first so we
#	perform the computation in a second run.
# # # # # # # # # # # # # # # # # # # # # # # # # # # 
import os

# change dir to the scripts directory
os.chdir( '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_v2' )

run_scripts = [ 'wrap_downscaler_cmip5_slurm.py',
				'wrap_downscaler_cmip5_slurm_tem.py',
				'wrap_downscaler_cru_40_slurm.py',
				'wrap_downscaler_cru_slurm.py',
				'wrap_downscaler_cru_slurm_tem.py',
				'wrap_downscaler_cru_ts40_slurm_tem.py' ]

for script in run_scripts:
	os.system( 'ipython {}'.format( script ) )

