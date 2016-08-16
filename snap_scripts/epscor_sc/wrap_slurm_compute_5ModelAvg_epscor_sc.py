# # # # # 
# wrap compute_5ModelAvg_epscor_sc.py for running on slurm
# # # # # 

def run_model( fn, base_dir, variable, model, scenario, units, metric ):
	head = '#!/bin/sh\n' + \
			'#SBATCH --ntasks=32\n' + \
			'#SBATCH --nodes=1\n' + \
			'#SBATCH --ntasks-per-node=32\n' + \
			'#SBATCH --account=snap\n' + \
			'#SBATCH --mail-type=FAIL\n' + \
			'#SBATCH --mail-user=malindgren@alaska.edu\n' + \
			'#SBATCH -p main\n'
	
	script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/compute_5ModelAvg_epscor_sc.py'
	with open( fn, 'w' ) as f:
		command = ' '.join([ 'ipython', script_path, '--', '-b', base_dir, '-v', variable ])
		f.writelines( head + "\n" + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1


if __name__ == '__main__':
	import os, subprocess
	
	variables = [ 'tas','pr','tasmin','tasmax' ]

	path = os.path.join( base_dir,'downscaled','slurm_log' )
	if not os.path.exists( path ):
		os.makedirs( path )

	os.chdir( path )
	for variable in variables:

		fn = os.path.join( path, 'slurm_run_compute_5ModelAvg_'+variable+'.slurm' )
		_ = run_model( fn, base_dir, variable, model, scenario, units, metric )
