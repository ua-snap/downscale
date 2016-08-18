# run all of them in one go with slurm

# # # EXAMPLE RUN SINGLE MODEL
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
	import os, subprocess

	script_directory = '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc'

	script_names = [ 'model_diffs_metrics_epscor_sc.py',
					'model_variability_metrics_epscor_sc_DECADAL_multidomain.py',
					'model_variability_metrics_epscor_sc_DECADAL_multidomain_cru.py' ]
	
	slurm_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/slurm_files'
	if not os.path.exists( slurm_path ):
		os.makedirs( slurm_path )
		
	os.chdir( slurm_path )

	for script in script_names:

		command = ' '.join(['ipython', os.path.join( script_directory, script ) ])

		fn = os.path.join( slurm_path, '_'.join(['slurm_run', script.replace('.py', '')]) + '.slurm' )
		_ = run_model( fn, command )