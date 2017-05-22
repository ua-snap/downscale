# # # # # 
# wrap compute_5ModelAvg_epscor_sc.py for running on slurm
# # # # # 

def run_model( fn, base_dir, variable, scenario ):
	head = '#!/bin/sh\n' + \
			'#SBATCH --ntasks=32\n' + \
			'#SBATCH --nodes=1\n' + \
			'#SBATCH --ntasks-per-node=32\n' + \
			'#SBATCH --account=snap\n' + \
			'#SBATCH --mail-type=FAIL\n' + \
			'#SBATCH --mail-user=malindgren@alaska.edu\n' + \
			'#SBATCH -p main\n'
	
	script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_v2/compute_5ModelAvg_cmip5.py'
	with open( fn, 'w' ) as f:
		command = ' '.join([ 'ipython', script_path, '--', '-b', base_dir, '-v', variable, '-s', scenario ])
		f.writelines( head + "\n" + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import os, subprocess
	
	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
	variables = [ 'pr', 'tasmin', 'tasmax', 'tas', 'hur', 'vap', 'clt' ]
	scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]

	slurm_path = os.path.join( base_dir, 'downscaled', 'slurm_log' )
	if not os.path.exists( slurm_path ):
		os.makedirs( slurm_path )
	
	os.chdir( slurm_path )

	for variable in variables:
		for scenario in scenarios:
			fn = os.path.join( slurm_path, 'slurm_run_compute_5ModelAvg_'+variable+'_'+scenario+'.slurm' )
			_ = run_model( fn, base_dir, variable, scenario )
