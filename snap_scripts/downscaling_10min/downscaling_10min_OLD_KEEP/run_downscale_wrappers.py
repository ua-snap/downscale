# run_downscale -- CMIP5 / CRU TS323 / CRU TS324_01 / CRU TS40
#  for variables: tas, pr 
#	perform the computation in a second run.
# # # # # # # # # # # # # # # # # # # # # # # # # # # 
import os

# change dir to the scripts directory
os.chdir( '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_10min' )

run_scripts = [ 'wrap_downscaler_cmip5_slurm.py',
				'wrap_downscaler_cru_40_slurm.py',
				'wrap_downscaler_cru_slurm.py',
				'wrap_downscaler_cru_324_01_slurm.py' ]

for script in run_scripts:
	os.system( 'ipython {}'.format( script ) )

