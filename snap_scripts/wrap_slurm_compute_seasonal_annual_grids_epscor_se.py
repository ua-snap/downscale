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
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_outputs_v2/decadal_monthly'
	ncpus = 32
	project = 'cmip5' # 'cru'
	variables = [ 'tasmin', 'tasmax', 'tas', 'pr' ]
	# variables = [ 'pr' ]
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4', '5ModelAvg' ]
	scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]

	for model in models:
		for scenario in scenarios:
			for variable in variables:
				# agg_metric = 'mean'
				if variable == 'pr':
					agg_metric = 'total'
				else:
					agg_metric = 'mean'
				
				if variable == 'tas' or variable == 'pr':
					# deal with different pathing to Matts downscaled variables on /Data
					if scenario == 'historical':
						base_path = '/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models'
					else:
						base_path = '/Data/Base_Data/Climate/AK_CAN_2km/historical/AR5_CMIP5_models'
				elif variable == 'tasmin' or variable == 'tasmax':
					base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_v2'
				
				# some stuff for output slurm log files
				path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_'+project+'_v2/slurm_log'
				
				if not os.path.exists( path ):
					os.makedirs( path )
				
				os.chdir( path )

				fn = os.path.join( path, 'slurm_run_compute_decadal_grids_'+'_'.join([variable, model, scenario, agg_metric])+'.slurm' )
				os.chdir( '/workspace/UA/malindgren/repos/downscale/snap_scripts' )
				command = ' '.join([ 'ipython', '/workspace/UA/malindgren/repos/downscale/snap_scripts/compute_decadal_grids_epscor_se.py', '--', '-b', base_path, '-o ', output_path, '-m ', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
				run_model( fn, command )

