# # # # # 
# wrap compute_5ModelAvg_epscor_sc.py for running on slurm
# # # # # 
def run_model( fn, model, scenario, variable ):
	head = '#!/bin/sh\n' + \
			'#SBATCH --ntasks=32\n' + \
			'#SBATCH --nodes=1\n' + \
			'#SBATCH --ntasks-per-node=32\n' + \
			'#SBATCH --account=snap\n' + \
			'#SBATCH --mail-type=FAIL\n' + \
			'#SBATCH --mail-user=malindgren@alaska.edu\n' + \
			'#SBATCH -p main,viz\n'
	
	script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_v2/format_for_alfresco.py'
	with open( fn, 'w' ) as f:
		command = ' '.join([ 'ipython', script_path, '--', '-m', model, '-s', scenario, '-v', variable ])
		f.writelines( head + "\n" + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import os, subprocess, itertools
	
	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
	models = [ 'GFDL-CM3', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'NCAR-CCSM4' ]
	scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	variables = ['tas', 'pr']

	slurm_path = os.path.join( base_dir, 'alfresco_1km', 'slurm_log' )
	if not os.path.exists( slurm_path ):
		os.makedirs( slurm_path )
	
	os.chdir( slurm_path )
	for model, scenario, variable in itertools.product( models,scenarios,variables ):
		fn = os.path.join( slurm_path, 'slurm_run_format_for_alfresco_{}_{}_{}.slurm'.format(model,scenario,variable) )
		_ = run_model( fn, model, scenario, variable )
	
	# CRU 
	models = [ 'CRU-TS40' ]
	scenarios = [ 'historical' ]

	slurm_path = os.path.join( base_dir, 'alfresco_1km', 'slurm_log' )
	if not os.path.exists( slurm_path ):
		os.makedirs( slurm_path )
	
	os.chdir( slurm_path )
	for model, scenario, variable in itertools.product( models,scenarios,variables ):
		fn = os.path.join( slurm_path, 'slurm_run_format_for_alfresco_{}_{}_{}.slurm'.format(model,scenario,variable) )
		_ = run_model( fn, model, scenario, variable )
	
	