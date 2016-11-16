# wrap downscaler for running on slurm
def run_model( fn, command ):
	import os, subprocess
	head = '#!/bin/sh\n' + \
			'#SBATCH --ntasks=32\n' + \
			'#SBATCH --nodes=1\n' + \
			'#SBATCH --ntasks-per-node=32\n' + \
			'#SBATCH --account=snap\n' + \
			'#SBATCH --mail-type=FAIL\n' + \
			'#SBATCH --mail-user=malindgren@alaska.edu\n' + \
			'#SBATCH -p main\n'
	
	with open( fn, 'w' ) as f:
		f.writelines( head + '\n' + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import subprocess, os
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/decadal_monthlies'
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4', '5ModelAvg' ]

	for model in models:
		
		# some stuff for output slurm log files
		path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/slurm_files_dot_dof/slurm_log'
		if not os.path.exists( path ):
			os.makedirs( path )
		
		os.chdir( path )

		fn = os.path.join( path, 'slurm_run_compute_decadal_grids_'+'_'.join([variable, model, scenario, agg_metric])+'.slurm' )
		os.chdir( '/workspace/UA/malindgren/repos/downscale/snap_scripts' )
		command = ' '.join([ 'ipython', '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/dot_dof_logs_cmip5_decadals.py', '--', 
							'-b', base_path, '-m ', model ])
		
		run_model( fn, command )

